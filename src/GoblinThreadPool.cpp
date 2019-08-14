#include "GoblinThreadPool.h"

namespace Goblin {

ThreadPool::ThreadPool(unsigned int coreNum,
    TLSManager* tlsManager):
    mTasksNum(0), mStartWork(false), mTLSManager(tlsManager) {
    mCoreNum = coreNum == 0 ?
        getMaxThreadNum() : std::min(getMaxThreadNum(), coreNum);
}


void ThreadPool::initWorkers() {
    if (mCoreNum == 1) {
        return;
    }
    for (size_t i = 0; i < mCoreNum; ++i) {
        mWorkers.push_back(
            new std::thread(&ThreadPool::taskEntry, this));
    }
}

void ThreadPool::taskEntry() {
    static thread_local TLSPtr tlsPtr;
    {
        std::unique_lock<std::mutex> lk(mStartMutex);
        while(!mStartWork) {
            mStartCondition.wait(lk);
        }
    }
    if (mTLSManager) {
        mTLSManager->initialize(tlsPtr);
    }
    while(true) {
        Task* task = nullptr;
        {
            std::lock_guard<std::mutex> lk(mTaskQueueMutex);
            if (mTasks.size() == 0) {
                break;
            }
            task= mTasks.back();
            mTasks.pop_back();
        }

        if (task) {
            task->run(tlsPtr);
        }

        {
            std::unique_lock<std::mutex> lk(mTaskQueueMutex);
            if (--mTasksNum == 0) {
                mTasksCondition.notify_all();
                break;
            }
        }
    }
    if (mTLSManager) {
        mTLSManager->finalize(tlsPtr);
    }
}

void ThreadPool::enqueue(const std::vector<Task*>& tasks) {
    if (mCoreNum == 1) {
        TLSPtr tlsPtr;
        mTLSManager->initialize(tlsPtr);
        for (size_t i = 0; i < tasks.size(); ++i) {
            tasks[i]->run(tlsPtr);
        }
        mTLSManager->finalize(tlsPtr);
        return;
    }

    if (mWorkers.size() == 0) {
        initWorkers();
    }
    
    {
        std::lock_guard<std::mutex> lk(mTaskQueueMutex);
        for (size_t i = 0; i < tasks.size(); ++i) {
            mTasks.push_back(tasks[i]);
        }
        mTasksNum += mTasks.size();
    }
};

void ThreadPool::waitForAll() {
    if (mCoreNum == 1) {
        return;
    }

    // let the worker start working
    {
        std::unique_lock<std::mutex> lk(mStartMutex);
        mStartWork = true;
    }
    mStartCondition.notify_all();

    // wake me up til all the taks finish
    { 
        std::unique_lock<std::mutex> lk(mTaskQueueMutex);
        while(mTasksNum != 0) {
            mTasksCondition.wait(lk);
        }
    } 
    cleanup();
}

void ThreadPool::cleanup() {
    for (size_t i = 0; i < mWorkers.size(); ++i) {
        if (mWorkers[i]->joinable()) {
            mWorkers[i]->join();
        }
        delete mWorkers[i];
    }
    mWorkers.clear();
    {
        std::unique_lock<std::mutex> lk(mStartMutex);
        mStartWork = false;
    }

}

}
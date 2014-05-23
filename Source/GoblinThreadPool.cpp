#include "GoblinThreadPool.h"

namespace Goblin {

    ThreadPool::ThreadPool(unsigned int coreNum):
        mTasksNum(0), mStartWork(false) {
        mCoreNum = coreNum == 0 ? 
            boost::thread::hardware_concurrency() : 
            min(boost::thread::hardware_concurrency(), coreNum);
    }


    void ThreadPool::initWorkers() {
        if(mCoreNum == 1) {
            return;
        }
        for(size_t i = 0; i < mCoreNum; ++i) {
            mWorkers.push_back(
                new boost::thread(&ThreadPool::taskEntry, this));
        }
    }

    void ThreadPool::taskEntry() {
        {
            boost::unique_lock<boost::mutex> lk(mStartMutex);
            while(!mStartWork) {
                mStartCondition.wait(lk);
            }
        }

        while(true) {
            Task* task = NULL;
            {
                boost::lock_guard<boost::mutex> lk(mTaskQueueMutex);
                if(mTasks.size() == 0) {
                    break;
                }
                task= mTasks.back();
                mTasks.pop_back();
            }

            if(task) {
                task->run();
            }

            {
                boost::unique_lock<boost::mutex> lk(mTaskQueueMutex);
                if(--mTasksNum == 0) {
                    mTasksCondition.notify_all();
                    break;
                }
            }

        }
    }

    void ThreadPool::enqueue(const vector<Task*>& tasks) {
        if(mCoreNum == 1) {
            for(size_t i = 0; i < tasks.size(); ++i) {
                tasks[i]->run();
            }
            return;
        }

        if(mWorkers.size() == 0) {
            initWorkers();
        }
    
        {
            boost::lock_guard<boost::mutex> lk(mTaskQueueMutex);
            for(size_t i = 0; i < tasks.size(); ++i) {
                mTasks.push_back(tasks[i]);
            }
            mTasksNum += mTasks.size();
        }
    };

    void ThreadPool::waitForAll() {
        if(mCoreNum == 1) {
            return;
        }

        // let the worker start working
        {
            boost::unique_lock<boost::mutex> lk(mStartMutex);
            mStartWork = true;
        }
        mStartCondition.notify_all();

        // wake me up til all the taks finish
        { 
            boost::unique_lock<boost::mutex> lk(mTaskQueueMutex);
            while(mTasksNum != 0) {
                mTasksCondition.wait(lk);
            }
        } 
        cleanup();
    }

    void ThreadPool::cleanup() {
        for(size_t i = 0; i < mWorkers.size(); ++i) {
            if(mWorkers[i]->joinable()) {
                mWorkers[i]->join();
            }
            delete mWorkers[i];
        }
        mWorkers.clear();
        {
            boost::unique_lock<boost::mutex> lk(mStartMutex);
            mStartWork = false;
        }

    }
}

#ifndef GOBLIN_THREAD_POOL_H
#define GOBLIN_THREAD_POOL_H

#include "GoblinThreadLocalStorage.h"
#include "GoblinUtils.h"
#include <thread>
#include <mutex>

namespace Goblin {

    class Task {
    public:
        virtual void run(TLSPtr& tls) = 0;
        virtual ~Task() {};
    };

    class ThreadPool {
    public:
        ThreadPool(unsigned int coreNum = 0,
            TLSManager* tlsManager = NULL);
        void enqueue(const vector<Task*>& tasks);
        void waitForAll();
        void cleanup();
    private:
        void initWorkers();
        void taskEntry();
    private:
        vector<std::thread*> mWorkers;
        unsigned int mCoreNum;

        std::condition_variable mTasksCondition;
        std::mutex mTaskQueueMutex;
        vector<Task*> mTasks;
        size_t mTasksNum;

        std::condition_variable mStartCondition;
        std::mutex mStartMutex;
        bool mStartWork;
        TLSManager* mTLSManager;
    };
}

#endif //GOBLIN_THREAD_POOL_H

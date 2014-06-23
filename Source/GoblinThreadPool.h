#ifndef GOBLIN_THREAD_POOL_H
#define GOBLIN_THREAD_POOL_H

#include <boost/thread.hpp>

#include "GoblinUtils.h"

namespace Goblin {

    class Task {
    public:
        virtual void run() = 0;
        virtual ~Task() {};
    };

    class ThreadPool {
    public:
        ThreadPool(unsigned int coreNum = 0);
        void enqueue(const vector<Task*>& tasks);
        void waitForAll();
        void cleanup();
    private:
        void initWorkers();
        void taskEntry();
    private:
        vector<boost::thread*> mWorkers;
        unsigned int mCoreNum;

        boost::condition_variable mTasksCondition;
        boost::mutex mTaskQueueMutex;
        vector<Task*> mTasks;
        size_t mTasksNum;

        boost::condition_variable mStartCondition;
        boost::mutex mStartMutex;
        bool mStartWork;
    };
}

#endif //GOBLIN_THREAD_POOL_H

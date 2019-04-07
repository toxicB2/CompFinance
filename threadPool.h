#pragma once

#include <future>
#include <thread>
#include "ConcurrentQueue.h"

using namespace std;

typedef packaged_task<bool(void)> Task;
typedef future<bool> TaskHandle;

class ThreadPool //singleton
{
	static ThreadPool theOnlyInstance;
    ConcurrentQueue<Task> workingQueue;
	vector<thread> allThreads;

    bool isActive;
	bool isInterrupted;

	static thread_local size_t protectedNumberOfGivenThread;

	void threadFunc(const size_t num)
	{
		protectedNumberOfGivenThread = num;
		Task t;

		while (!isInterrupted) 
		{
			workingQueue.pop(t);
			if (!isInterrupted) t();			
		}
	}

    ThreadPool() : isActive(false), isInterrupted(false) {}

public:
	//Getters
	static ThreadPool* getInstance() { return &theOnlyInstance; }
	size_t getTotalNumerOfThreads() const { return allThreads.size(); }
	static size_t getThreadNumer() { return protectedNumberOfGivenThread; }

	void start(const size_t nThread = thread::hardware_concurrency() - 1)
	{
        if (!isActive)  //  Only start once
        {
            allThreads.reserve(nThread);

            //	Launch threads on threadFunc and keep handles in a vector
            for (size_t i = 0; i < nThread; i++)
                allThreads.push_back(thread(&ThreadPool::threadFunc, this, i + 1));

            isActive = true;
        }
	}

	//	Destructor
    ~ThreadPool()
    {
        stop();
    }
        
    void stop()
	{
        if (isActive)
        {
            //	Interrupt mode
            isInterrupted = true;

            //	Interrupt all waiting threads
            workingQueue.interrupt();

            //	Wait for them all to join
            for_each(allThreads.begin(), allThreads.end(), mem_fn(&thread::join));

            //  Clear all threads
            allThreads.clear();

            //  Clear the queue and reset interrupt
            workingQueue.clear();
            workingQueue.resetInterrupt();

            //  Mark as inactive
            isActive = false;

            //  Reset interrupt
            isInterrupted = false;
        }
	}

	//	Forbid copies etc
	ThreadPool(const ThreadPool& rhs) = delete;
	ThreadPool& operator=(const ThreadPool& rhs) = delete;
	ThreadPool(ThreadPool&& rhs) = delete;
	ThreadPool& operator=(ThreadPool&& rhs) = delete;

	//	Spawn task
	template<typename Callable>
	TaskHandle spawnTask(Callable c)
	{
		Task t(move(c));
		TaskHandle f = t.get_future();
		workingQueue.push(move(t));
		return f;
	}

	//	Run queued tasks synchronously 
	//	while waiting on a future, 
	//	return true if at least one task was run
	bool activeWait(const TaskHandle& f)
	{
		Task t;
		bool b = false;

		//	Check if the future is ready without blocking
		//	The only syntax C++11 provides for that is
		//	wait 0 seconds and return status
		while (f.wait_for(0s) != future_status::ready)
		{
			//	Non blocking
			if (workingQueue.tryPop(t)) 
			{
				t();
				b = true;
			}
			else //	Nothing in the queue: go to sleep
			{
				f.wait();
			}
		}

		return b;
	}
};
#pragma once

#include <queue>
#include <mutex>
using namespace std;

template <class T>
class ConcurrentQueue
{
    queue<T> workingQueue;
	mutable mutex mutexObj;
	condition_variable condVar;
	bool isInterrupted;

public:
	ConcurrentQueue() : isInterrupted(false) {}
	~ConcurrentQueue() { interrupt(); }

	bool empty() const
	{
		lock_guard<mutex> lk(mutexObj);
		return workingQueue.empty();
	}

	bool tryPop(T& t)
	{
		lock_guard<mutex> lk(mutexObj);
		if (workingQueue.empty()) return false;
		t = move(workingQueue.front());
		workingQueue.pop();

		return true;
	}

	void push(T t) //	Pass t byVal or move with push( move( t))
	{
		{
			lock_guard<mutex> lk(mutexObj);
			workingQueue.push(move(t));
		}	
        //	Unlock before notification 
		condVar.notify_one();
	}

	//	Wait if empty
	bool pop(T& t)
	{
		unique_lock<mutex> lk(mutexObj);
		//	Wait if empty, release lock until notified 
		while (!isInterrupted && workingQueue.empty()) condVar.wait(lk);

		if (isInterrupted) 
			return false;

		t = move(workingQueue.front());
		workingQueue.pop();

		return true;

	}	

	void interrupt()
	{
        {
            lock_guard<mutex> lk(mutexObj);
            isInterrupted = true;
        }
		condVar.notify_all();
	}

    void resetInterrupt(){ isInterrupted = false;}

    void clear()
    {
        queue<T> empty;
        swap(workingQueue, empty);
    }
};
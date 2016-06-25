/*******************************************************************************
 * gpuATTRACT framework
 * Copyright (C) 2015 Uwe Ehmann
 *
 * This file is part of the gpuATTRACT framework.
 *
 * The gpuATTRACT framework is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The gpuATTRACT framework is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

#ifndef SAVEWORKQUEUE_H_
#define SAVEWORKQUEUE_H_

#include <queue>
#include <mutex>
#include <condition_variable>
#include <atomic>

namespace as {


template <class T>
class ThreadSafeQueue {
public:
	ThreadSafeQueue() :
		_terminate(false)
	{}

	~ThreadSafeQueue() {}

	T waitAndPop () {
		std::unique_lock < std::mutex > uniqueLock(_mutex);
		while( (_queue.size() == 0) ) {
			if (_terminate.load()) {
				return T();
			}
			_condVar.wait(uniqueLock);
		}
		T item = _queue.front();
		_queue.pop();
		return item;
	}

	bool tryPop(T& value) {
		std::lock_guard<std::mutex> guard(_mutex);
		if (_queue.size() == 0) {
			return false;
		}
		value = _queue.front();
		_queue.pop();
		return true;
	}

	void push(T const& value) {
		std::lock_guard<std::mutex> guard(_mutex);
		_queue.push(value);
		_condVar.notify_one();
	}

	void signalTerminate() {
		_terminate = true;
		_condVar.notify_one();
	}

	bool terminates() const {
		return _terminate.load();
	}

	unsigned size() const {
		std::lock_guard<std::mutex> guard(_mutex);
		return _queue.size();
	}

	bool empty() const {
		std::lock_guard<std::mutex> guard(_mutex);
		return _queue.empty();
	}

private:
	std::queue<T> _queue;
	mutable std::mutex _mutex;
	std::condition_variable _condVar;

	std::atomic<bool> _terminate;
};

} // namespace

#endif /* SAVEWORKQUEUE_H_ */

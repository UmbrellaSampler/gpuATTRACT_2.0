/*
 * RequestStatus.h
 *
 *  Created on: Feb 8, 2016
 *      Author: uwe
 */

#ifndef SRC_REQUESTSTATUS_H_
#define SRC_REQUESTSTATUS_H_

namespace as {

class RequestStatus {
private:
	enum Status {
		InProgress,
		Completed,
		invalid
	};

public:

	RequestStatus() :
		_status(invalid) {}

	void setInProgress() noexcept {
		_status = InProgress;
	}

	void setCompleted() noexcept {
		_status = Completed;
	}

	bool inProgress() noexcept {
		return _status == InProgress;
	}

	bool completed() noexcept {
		return _status == Completed;
	}

private:
	Status _status;
};

}  // namespace as



#endif /* SRC_REQUESTSTATUS_H_ */

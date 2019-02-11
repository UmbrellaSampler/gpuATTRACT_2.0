/*
 * Test_Server.cpp
 *
 *  Created on: Dec 30, 2015
 *      Author: uwe
 */

#include <gtest/gtest.h>
#include <numeric>

#include "Server.h"
#include "Request.h"
#include "Allocator.h"

#include "Service_Mock.h"
#include "Service_TimeOut.h"

#include "HostAllocator.tpp"
#include <memory>
//#include <iostream>
//using std::cout;
//using std::endl;

using ::testing::_;
using ::testing::AtLeast;
using ::testing::Invoke;

namespace test {

TEST(Server, Construction) {

	try {
		as::Server<Service_Mock::genericTypes_t> server(nullptr);
		FAIL();
	} catch (std::exception& e) {
		ASSERT_STREQ( "Invalid service (nullptr).", e.what() );
	}

	auto service = std::make_shared<Service_Mock>();

	try {
		as::Server<Service_Mock::genericTypes_t> server(service);
		FAIL();
	} catch (std::exception& e) {
		ASSERT_STREQ( "Invalid service allocators (nullptr).", e.what() );
	}

	try {
		as::Server<Service_Mock::genericTypes_t> server(service);
		EXPECT_TRUE(true);
	} catch (std::exception& e) {
		FAIL();
	}
}

TEST(Server, Interface_setService) {

	auto service = std::make_shared<Service_Mock>();
	as::Server<Service_Mock::genericTypes_t> server(service);
	try {
		server.setService(nullptr);
		FAIL();
	} catch (std::exception& e) {
		ASSERT_STREQ( "Invalid service (nullptr).", e.what() );
	}

	try {
		server.setService(service);
		FAIL();
	} catch (std::exception& e) {
		ASSERT_STREQ( "Invalid service allocators (nullptr).", e.what() );
	}

	try {
		server.setService(service);
		EXPECT_TRUE(true);
	} catch (std::exception& e) {
		FAIL();
	}
}

TEST(Server, Interface_submit) {
	using request_t = as::Server<Service_Mock::genericTypes_t>::request_t;
	auto serviceMock = std::make_shared<Service_Mock>();
	as::Server<Service_Mock::genericTypes_t> server(serviceMock);
	ON_CALL(*serviceMock, createItemProcessor())
		.WillByDefault(Invoke(serviceMock.get(), &Service_Mock::createItemProcessorImpl));
	EXPECT_CALL(*serviceMock, createItemProcessor())
		.Times(AtLeast(0));
	server.createWorkers(1);

	auto request = std::make_shared<request_t>();
	try {
		server.submit(request);
		FAIL();
	} catch (std::exception& e) {
		ASSERT_STREQ( "Invalid input buffer (nullptr).", e.what() );
	}

	request.reset();
	int size = 100000;
	int input[size];
	int result[size];
	request->setDataAndSize(input, size);
	try {
		server.submit(request);
		request->setDataAndSize(input, size);
		server.submit(request);
		FAIL();
	} catch (std::exception& e) {
		ASSERT_STREQ( "Request has already been submitted.", e.what() );
		server.wait(request, result); // to ensure that buffer return to buffer manager
	}

	request.reset();
	try {
		for (int i = 0; i < 10; ++i) {
			request->setDataAndSize(input, 0);
			server.submit(request);
			EXPECT_TRUE(true);
			server.wait(request, result);
		}
	} catch (std::exception& e) {
		FAIL();
	}

}

TEST(Server, Interface_wait) {
	using request_t = as::Server<Service_Mock::genericTypes_t>::request_t;
	auto serviceMock = std::make_shared<Service_Mock>();
	as::Server<Service_Mock::genericTypes_t> server(serviceMock);
	ON_CALL(*serviceMock, createItemProcessor())
		.WillByDefault(Invoke(serviceMock.get(), &Service_Mock::createItemProcessorImpl));
	EXPECT_CALL(*serviceMock, createItemProcessor())
		.Times(AtLeast(0));
	server.createWorkers(1);
	auto request = std::make_shared<request_t>();
	request->reset();
	int size = 1;
	int input[size];
	int result[size];

	request->setDataAndSize(input, size);
	try {
		server.submit(request);
		auto otherRequest = std::make_shared<request_t>();
		server.wait(otherRequest, result);
		FAIL();
	} catch (std::exception& e) {
		ASSERT_STREQ( "Either request has not been submitted or clientBuffer invalid (nullptr).", e.what() );
		server.wait(request, result); // to ensure that buffer return to buffer manager
	}

	request->setDataAndSize(input, size);
	try {
		server.submit(request);
		server.wait(request, nullptr);
		FAIL();
	} catch (std::exception& e) {
		ASSERT_STREQ( "Either request has not been submitted or clientBuffer invalid (nullptr).", e.what() );
		server.wait(request, result); // to ensure that buffer return to buffer manager
	}

	using request_t = as::Server<Service_TimeOut>::request_t;

	as::Server<Service_TimeOut::genericTypes_t> server_to(std::make_shared<Service_TimeOut>());
	server_to.createWorkers(1);
	server_to.setWaitTime(10);

	auto request_to = std::make_shared<request_t>();
	request_to->setDataAndSize(input, size);
	try {
		server_to.submit(request_to);
		server_to.wait(request_to, result);
		FAIL();
	} catch (std::exception& e) {
		ASSERT_STREQ( "Request takes too long to finish for unknown reasons", e.what() );
	}

}

TEST(Server, Interface_setItemSize) {
	auto serviceMock = std::make_shared<Service_Mock>();
	as::Server<Service_Mock::genericTypes_t> server(serviceMock);

	try {
		server.setItemSize(0);
		FAIL();
	} catch (std::exception& e) {
		ASSERT_STREQ( "Item size must be greater than zero.", e.what() );
	}

	try {
		server.setItemSize(1368);
		EXPECT_TRUE(true);
	} catch (std::exception& e) {
		FAIL();
	}
}


constexpr unsigned maxTestSize = 100000;
constexpr unsigned testItemSize = 1000;
constexpr unsigned testCommon = 2;

class ServerIntegration: public ::testing::Test {
protected: // could also be public according to gtest Primer.md
	ServerIntegration() :
		sh_serviceMock(std::make_shared<Service_Mock>()),
		server(sh_serviceMock),
		common(testCommon){}

	virtual void SetUp() {
		initInput();
	}

	virtual void TearDown() {}

	using server_t = as::Server<Service_Mock::genericTypes_t>;
	using request_t = as::Server<Service_Mock::genericTypes_t>::request_t;

	request_t createRequest(unsigned inputSize) {
		if (inputSize > maxTestSize) {
			inputSize = maxTestSize;
		}
		request_t request;
		request.setDataAndSize(input, inputSize);
		request.setCommon(common);
		return request;
	}

	std::shared_ptr<Service_Mock> sh_serviceMock; // mock service to be injected into the service
	server_t server;				// associated server

	int common;			// the common element of the request

	Service_Mock::input_t input[maxTestSize];
	Service_Mock::result_t result[maxTestSize];

private:

	void initInput() {
		std::iota(input, input + maxTestSize, 0);
	}
};


TEST_F(ServerIntegration, WithServiceMock)
{
	ON_CALL(*sh_serviceMock, createItemProcessor())
		.WillByDefault(Invoke(sh_serviceMock.get(), &Service_Mock::createItemProcessorImpl));
	EXPECT_CALL(*sh_serviceMock, createItemProcessor())
		.Times(AtLeast(1));

	std::vector<unsigned> itemSizes = {1,2,3,5,7,55,111,390,557,1000};
	std::vector<unsigned> inputSizes = {10,537,1001,1,2,3, 100000};
	std::vector<unsigned> poolSize = {1,2,3,4};

	for (unsigned workers : poolSize) {
		server.createWorkers(workers);

		for (unsigned inputSize : inputSizes) {
			for (unsigned itemSize : itemSizes) {
				auto request  = std::make_shared<request_t>(createRequest(inputSize));
				server.setItemSize(itemSize);
				server.submit(request);
				server.wait(request, result);
				for(unsigned i = 0; i < inputSize; ++i) {
					EXPECT_EQ(result[i], input[i]*common);
				}
			}
		}
	}

}

} // namespace test














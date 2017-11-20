/*
 * CompareData.h
 *
 *  Created on: Nov 12, 2017
 *      Author: glenn
 */

#ifndef COMPAREDATA_H_
#define COMPAREDATA_H_


#include <iostream>
#include <string>

template <typename REAL>
struct evalData{
	int size;
	int *index;
	REAL *data;
	evalData(int dataSize) : index(new int[dataSize]), data(new REAL[dataSize]),size(dataSize) {}
};

template <typename REAL>
class CompareData{
private:
	REAL* m_referenceData;
	REAL* m_testData;//=(REAL*) malloc(875*sizeof(REAL));;
	int m_dataSize;
	REAL m_epsilon;
	evalData<REAL> m_result;
//int dataSize,
public:
	 //explicit MyArray(std::size_t s) : buf(new float[s]), size(s) { }
	explicit CompareData(const int s, REAL epsilon): m_dataSize(s),m_result(s),m_referenceData(new REAL[s]),m_testData(new REAL[s])
	{
		//:m_dataSize(dataSize){
		//m_dataSize=s;
		 //m_referenceData=new REAL[m_dataSize];
		 //m_testData=
//		 m_result.index=new int[m_dataSize];
//		 m_result.data=new REAL[m_dataSize];
		m_epsilon=epsilon;

//		for(int i=0;i<m_dataSize;i++){
//			m_result.index[i]=(float)i;
//						//std::cout <<i<<std::endl;
//						}
//
//		for(int i=0;i<m_dataSize;i++){
//				std::cout <<m_result.index[i]<<endl;
//				}
	}
	~CompareData(){
			delete [] m_referenceData;
			delete [] m_testData;
			delete [] m_result.index;
			delete [] m_result.data;
	}

	void setEpsilon(REAL epsilon){
		m_epsilon=epsilon;
	}

	REAL getEpsilon(){
			return m_epsilon;
	}
	REAL* referenceData(){
		return m_referenceData;
	}
	REAL* testData(){
		return m_testData;
	}



	void evaluateRatio(){
		int numFailure=0;
		float ratio;
		for(int i=0;i<m_dataSize;i++){
			ratio=m_testData[i]/m_referenceData[i];
			if (!(1-m_epsilon < ratio && 1+m_epsilon >ratio)){	m_result.data[numFailure]=ratio;	m_result.index[numFailure]=i;  numFailure++;}
		}
		m_result.size=numFailure;
	}


	void evaluateDifference(){
			int numFailure=0;
			float delta;
			for(int i=0;i<m_dataSize;i++){
				delta=m_testData[i]-m_referenceData[i];
				if (!(delta*delta < m_epsilon)){	m_result.data[numFailure]=delta;	m_result.index[numFailure]=i;  numFailure++;}
			}
			m_result.size=numFailure;
		}

	void evaluateExact(){
		int numFailure=0;
		float ratio;
		for(int i=0;i<m_dataSize;i++){
			ratio=m_testData[i]/m_referenceData[i];
			if ( !(ratio == 1.0)){	m_result.data[numFailure]=ratio;	m_result.index[numFailure]=i;	numFailure++;}
		}
		m_result.size=numFailure;
		}

	float getAccuracy(){
		return ((float)m_dataSize-(float)m_result.size)/(float)m_dataSize;
	}

	float getNumFailure(){
		return m_result.size;
	}
	evalData<REAL>* getResult(){
		return *m_result;
	}

	void writeResultToFile(std::string filename){
		using namespace std;

		ofstream myfile;
		myfile.open (filename);
		myfile << "## numFailure" << endl;
		myfile << "#" << m_result.size << endl;

		myfile << "## Accuracy" << endl;
		myfile << "#" << ((float)m_dataSize-(float)m_result.size)/(float)m_dataSize<<endl;

		myfile << "## epsilon" << endl;
		myfile << "#" << m_epsilon <<endl;

		myfile << "## result" << endl;
		myfile << "### indices \t ratio(testdata/referencedata\t testData\t referenceData"<< endl;

		for(int i=0;i<m_result.size;i++){
			myfile << m_result.index[i] << "\t" << m_result.data[i]<< "  " <<  m_testData[m_result.index[i]]<< "   " << m_referenceData[m_result.index[i]] << endl;
		}
		myfile.close();

		}

};


#endif /* COMPAREDATA_H_ */

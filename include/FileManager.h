/*
 * FileManager.h
 *
 *  Created on: Dec 27, 2011
 *      Author: user
 */
#include <iostream>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include "Constant.h"
#ifndef FILEMANAGER_H_
#define FILEMANAGER_H_
using namespace std;
class FileManager {

	ifstream f;
	ofstream of;
public:
	FileManager();
	FileManager(const char *,Constant::MODE m);
	virtual ~FileManager();

	void readFile(double *result);
	void readFile(std::vector<int>& result);
	void writeFile(const double t);
	void writeFile(const double *w,int num);
	void close();
};

#endif /* FILEMANAGER_H_ */

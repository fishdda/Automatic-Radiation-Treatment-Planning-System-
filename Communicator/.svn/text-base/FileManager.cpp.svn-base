/*
 * FileManager.cpp
 *
 *  Created on: Dec 27, 2011
 *      Author: user
 */

#include "FileManager.h"

FileManager::FileManager()
{

}
FileManager::FileManager(const char* fileName, Constant::MODE m) {
	// TODO Auto-generated constructor stub
	//of.open(fileName);
	if (m == Constant::Read) {
		f.open(fileName);
		if (!f.is_open())
			printf("Can't open file\n");
	} else if (m == Constant::Write) {
		of.open(fileName);
		if (!of.is_open())
			printf("Can't open file\n");
	}

}
FileManager::~FileManager() {

}
void FileManager::close() {
	f.close();
	of.close();
}
void FileManager::readFile(double *result) {
	char line[50];
	int k = 0;

	while (!f.eof()) {
		//cout <<line <<endl;
		f >> line;
		//cout<<line<<endl;
		result[k] = atof(line);
		k++;
	}

}
void FileManager::readFile(std::vector<int>& result) {
	char line[50];
	int k = 0;

	while (!f.eof()) {
		//cout <<line <<endl;
		f >> line;
		// cout<<line<<endl;
		//cout<<line<<endl;
		int id = atoi(line);
		result.push_back(id);

	}

}
void FileManager::writeFile(const double t) {

	of << t << "\n";

}
void FileManager::writeFile(const double *w, int num) {

	for (int i = 0; i < num; i++) {
		of << w[i] << "\n";
	}
}

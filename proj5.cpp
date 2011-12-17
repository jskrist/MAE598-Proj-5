#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <math.h>
#include <mpi.h>
#include <string>
#include "parser.h"

using namespace std;

const int MAXCHARS = 80;
std::string szComment = "*~*";

double update(int &x,int &y,int &z);
void   setDims(int &nSize);
void   getLineNums(const MPI::File dataFile, int nlineNum[], int maxLine);
int	   getNumLines(const MPI::File dataFile);

//**************************************************
//  Global Variables
//**************************************************

CParser Parser;

//**************************************************
//  Main Program
//**************************************************
int main(int argc, char** argv)
{
	//initialize MPI
	MPI::Init();

	// Rank and size of processor layout
	int rank = 0;
	int size = 0;

	// Number of lines in the input file
	int numLines = 0;

	//Start and stop points in file
	int start = 0, stop = 0;

	//Vars for word stats
	string buf;
	int wordCnt = 0, nTokens = 0;
	vector<string> wordVect;
	vector<string> tmpWordVect;
	vector<string>::iterator it;
	it = wordVect.begin();

	vector<int>    cntVect;
	vector<double> probVect;

	//Parser
    CParser Parse;

	//Line number
	int nLineNum = 0;

	//MPI Communicator
	Comm = MPI_COMM_WORLD;

	rank = MPI::COMM_WORLD.Get_rank();
    size = MPI::COMM_WORLD.Get_size();

	//Define MPI String type
	MPI_Datatype MPI_STR;
	MPI_Type_contiguous(16,MPI::CHAR,&MPI_STR);
	MPI_Type_commit(&MPI_STR);

	// Open file to which data will be written
	MPI::File dataFile;
	dataFile = MPI::File::Open( oldComm, 
								"/scratch/jskristo/corpus.txt",
								MPI::MODE_CREATE | MPI::MODE_RDWR,
								MPI::INFO_NULL);

	int fileSize = dataFile.Get_size();

	if (fileSize == 0)
	{
		cout << "File is empty no stats available.\n";
		return 1;
	}
	else
	{
		if(rank == 0)
			numLines = getNumLines(dataFile);

		MPI_Allreduce( MPI_IN_PLACE, &numLines, 1, MPI_INT, MPI_MAX, Comm);

		int* nlineNum;
		nLineNum = (int*) malloc (numLines);
		if (nLineNum == NULL)
			exit (1);

		if(rank == 0)
			getLineNums(dataFile, nlineNum, numLines);

		MPI_Allreduce( MPI_IN_PLACE, &nLineNum, numLines, MPI_INT, MPI_MAX, Comm);

		if(rank <= numLines)
		{
			start = floor(numLines / size) * rank;
			if(rank = size - 1)
				stop = numLines;
			else
				stop = ( floor(numLines / size) * (rank + 1) ) - 1;
		}
	}

	dataFile.Seek(nlineNum[start], MPI_SEEK_SET);

	szTmpStr = "";
	for (int j = start; j < stop; j++)
	{
		if (!Parse.GetTokens (dataFile, tmpWordVect, nTokens, " "))
		{
			cout << "error getting tokens\n";
			return 1;
		}

		for(int k = 0; k < nTokens; k++)
		{
			if(wordVect.size() == 0)
			{
				wordVect.insert(it, tmpWordVect[k]);
				cntVect.insert(it, 1);
			}
			else
			{
				for(int l = 0; l < wordVect.size(); l++)
				{
					if (tmpWordVect[k] < wordVect[l]);
					{
						if(l == wordVect.size() - 1)
						{
							wordVect.insert(it + l, tmpWordVect[k]);
							cntVect.insert(it + l, 1);
							break;
						}
						else if( tmpWordVect[k] > wordVect[l + 1] )
						{
							wordVect.insert(it + l, tmpWordVect[k]);
							cntVect.insert(it + l, 1);
							break;
						}
					}
					else if (tmpWordVect[k] == wordVect[l])
					{
						cntVect[l]++;
						break;
					}
					else
					{
						wordVect.push_back(tmpWordVect[k]);
						cntVect.push_back(1);
						break;
					}
				}
			}
		}

		dataFile.Iread(&buf, wordCnt, MPI_STR);
	}

	dataFile.Close();

	MPI::Finalize();
	return 0;
}

void   getLineNums(const MPI::File dataFile, int vnlineNum[], int maxLine)
{
	string szTmpStr = "";
    CParser Parse;
	int nLineNum = 0;

    // read the file once
	for(int i = 0; i < maxLine; i++)
	{
		if (!Parse.ReadNextLine (dataFile, nLineNum, szTmpStr, MAXCHARS, szComment, true))
			cout << "Error Parsing line " << nLineNum << "\n";

		vnLineNum[nLineNum - 1] = tellg();
	}

	Rewind(dataFile);
}

int	   getNumLines(const MPI::File dataFile)
{
	string szTmpStr = "";
    CParser Parse;
	int nLineNum = 0;

    // read the file once
    if (!Parse.ReadNextLine (dataFile, nLineNum, szTmpStr, MAXCHARS, szComment, true))
		cout << "Error Parsing line " << nLineNum << "\n";

	//while not at end of file, continue reading file
	while (!dataFile.eof( ))
	{
		if (!Parse.ReadNextLine (dataFile, nLineNum, szTmpStr, MAXCHARS, szComment, true))
			cout << "Error Parsing line " << nLineNum << "\n";
	}
	Rewind(dataFile);
	return nLineNum;
}

void   setDims(int &nSize)
{
	int tmpDim = (int)floor(pow((double)nSize,0.33333333333333333333));
	int tmpSize = nSize / tmpDim;

	dims[0] = tmpDim;

	tmpDim = (int)floor(pow((double)tmpSize,0.5));

	dims[1] = tmpDim;

	tmpDim = tmpSize / tmpDim;

	dims[2] = tmpDim;	
}

double update(int &x, int &y, int &z)
{
	double result = 0;

	double xP1=0,
		   xM1=0,
		   yP1=0,
		   yM1=0,
		   zP1=0,
		   zM1=0,
		   nCur=nodes[z][y][x];

	if(x+1 < i)
		xP1=nodes[z][y][x+1];
	else
		xP1=nCur;

	if(x-1 >= 0)
		xM1=nodes[z][y][x-1];
	else
		xM1=nCur;

	if(y+1 < j)
		yP1=nodes[z][y+1][x];
	else
		yP1=nCur;

	if(y-1 >= 0)
		yM1=nodes[z][y-1][x];
	else
		yM1=nCur;

	zP1=nodes[z+1][y][x];
	zM1=nodes[z-1][y][x];

	result = alpha*dt*(pow((double)dw,2)*(xM1-2*nCur+xP1)
					 + pow((double)dh,2)*(yM1-2*nCur+yP1)
					 + pow((double)dl,2)*(zM1-2*nCur+zP1)) + nCur;

	return result;
}




#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include<iterator>
#include <algorithm>
#include<mpi.h>
#include <iomanip>

using namespace std;
const double KM = 1000.0;
const int TAG_DATA = 0, TAG_QUIT = 1;
const string IN_FILE_FOODBANKS = "foodbanks.dat";
const string IN_FILE_RESIDENCES = "residences.dat";

typedef struct Counter
{
	int range1,range2,range3,range4;
	double preset1 ,preset2 ,preset3 ,preset4;
}  countType;

MPI_Datatype createRecType()
{
	// Set-up the arguments for the type constructor
	MPI_Datatype new_type;
	int count = 2;

	int blocklens[] = { 4, 4 };
	MPI_Aint indices[2];
	indices[0] = 0;
	MPI_Type_extent( MPI_INT, &indices[1] );
	indices[1] *= 4;    // There are 2 doubles
	MPI_Datatype old_types[] = { MPI_INT ,MPI_DOUBLE};
	// Call the data type constructor
	MPI_Type_struct(count, blocklens, indices, old_types, &new_type);
	MPI_Type_commit(&new_type);

	return new_type;
}
struct Foodbank
{
	double x;
	double y;


	//input stream overload
	friend std::istream &operator>>(std::istream &is, Foodbank &f) {
		return is>>f.x>>f.y;
	}
	//output stream overload
	friend std::ostream &operator<<(std::ostream &os, Foodbank const &f) {
		return os << f.x  << "\t"<< f.y<<"\t";
	}
};
struct Residence
{
	double  x;
	double y;
	//input stream overload
	friend std::istream &operator>>(std::istream &is, Residence &r) {
		return is>>r.x>>r.y;
	}
	//output stream overload
	friend std::ostream &operator<<(std::ostream &os, Residence const &r) {
		return os << r.x  << "\t"<< r.y<<"\t";

	}
};


double distanceCalculate(double x1, double y1, double x2, double y2)
{
	double x = x1 - x2;
	double y = y1 - y2;
	double dist;

	dist = pow(x,2)+pow(y,2);           //calculating distance by euclidean formula
	dist = sqrt(dist);                  //sqrt is function in math.h

	return dist;
}
//get short destainsce
double getShortestDistances(Residence rs,vector<Foodbank> fs)
{
	double distances;



	std::vector<double>tempVec;
	//for(auto f:fs)

	std::vector<Foodbank>::iterator iter;

	for( iter = fs.begin(); iter != fs.end(); iter++ )
	{
	

		tempVec.push_back( distanceCalculate(iter->x,iter->y,rs.x,rs.y) / KM);

	}
	distances = *std::min_element(tempVec.begin(), tempVec.end());


	return distances;
}


void analysis_range(double d,countType &ctp)
{

	//from 0 to 1.00 KM
	if(d >= 0.00 && d <= 1.00)
	{
		ctp.range1 +=1;
		ctp.preset1 = ((double) ctp.range1/(ctp.range1+ctp.range2+ctp.range3+ctp.range4)) * 100.00;
	}
	//from 1.00 to 2.00 KM
	else if(d > 1.00 && d<= 2.00)
	{
		ctp.range2 +=1;
		ctp.preset2 =  ((double) ctp.range2/(ctp.range1+ctp.range2+ctp.range3+ctp.range4)) * 100.00;
	}
	//from 1.00 to 2.00 KM
	else if (d > 2.00 && d<= 5.00)
	{
		ctp.range3 +=1;
		ctp.preset3 = ((double) ctp.range3/(ctp.range1+ctp.range2+ctp.range3+ctp.range4)) * 100.00;
	}
	//grater than 5.00 KM
	else if (d > 5.00 )
	{
		ctp.range4 +=1;
		ctp.preset4 = ((double) ctp.range4/(ctp.range1+ctp.range2+ctp.range3+ctp.range4)) * 100.00;
	}

}


void processData(int rank,int numProcs)
{

	static countType count;
	MPI_Datatype recType = createRecType();

	//read file and populate the vectors
	ifstream foodbankFile("foodbanks.dat");
	ifstream residenceFile("residences.dat");

	// populate datavector
	std::vector<Foodbank> foodbankData((std::istream_iterator<Foodbank>(foodbankFile)),
		std::istream_iterator<Foodbank>());

	Residence res;
	int numLines = 0;


	while(!residenceFile.eof())
	{
		residenceFile >> res.x >>res.y;


		if ( numLines % numProcs == rank)
		{

			analysis_range(getShortestDistances(res,foodbankData),count);

		}
		++numLines;

	}

	countType* countBuff  = new countType[numProcs];
	MPI_Gather(&count, 1, recType, countBuff, 1, recType,0, MPI_COMM_WORLD);

	/*std::cout<< "for Rank"<<rank<< ",from 0 to 1.00 KM:"<<count.range1<<",%"<<count.preset1
	<<",from 1.00 to 2.00 KM:"<<count.range2<<",%"<<count.preset2<<",from 2.00 to 5.00 KM:"
	<<count.range3<<",%"<<count.preset3<<",grater than 5.00 KM:"<<count.range4<<",%"<<count.preset3<<std::endl;*/
	if(rank == 0)
	{

		static countType countArggResult;
		for (int p = 0; p < numProcs; ++p)
		{

			countArggResult.range1 +=countBuff[p].range1;
			countArggResult.range2 +=countBuff[p].range2;
			countArggResult.range3 +=countBuff[p].range3;
			countArggResult.range4 +=countBuff[p].range4;

			int proseRank  = p+1;
			int totalNumberOFAddress = countBuff[p].range1 + countBuff[p].range2 +countBuff[p].range3 + countBuff[p].range4;
			cout<<"For Rank"<<proseRank<<endl;
			cout<<"Total Number OF Address ="<<totalNumberOFAddress<<endl;
		cout<<"Nerest Foodbank\t"<<"# of Address"<<"\t %Address"<<endl;
		cout<<"---------------\t"<<"------------"<<"\t---------"<<endl;
			cout<<"0.00 - 1.00 \t"<<countBuff[p].range1<<"\t\t%"<<countBuff[p].preset1<<endl;
			cout<<"1.00 - 2.00 \t"<<countBuff[p].range2<<"\t\t %"<<countBuff[p].preset2<<endl;
			cout<<"2.00 - 5.00 \t"<<countBuff[p].range3<<"\t\t %"<<countBuff[p].preset3<<endl;
			cout<<"    5.00 \t"<<countBuff[p].range4<<"\t\t %"<<countBuff[p].preset4<<endl;
			cout<<endl;
			cout<<endl;

		}
		cout<<endl;
		cout<<endl;
		// calulate the total 
		int totalNumberOFAddressArrg = countArggResult.range1 + countArggResult.range2+ countArggResult.range3+ countArggResult.range4;
		countArggResult.preset1 = ((double)countArggResult.range1/ totalNumberOFAddressArrg) * 100.00;
		countArggResult.preset2 = ((double)countArggResult.range2/ totalNumberOFAddressArrg) * 100.00;
		countArggResult.preset3 = ((double)countArggResult.range3/ totalNumberOFAddressArrg) * 100.00;
		countArggResult.preset4 = ((double)countArggResult.range4/ totalNumberOFAddressArrg) * 100.00;
		cout<<"total Number OF Address Aggregation ="<<totalNumberOFAddressArrg<<endl;
		cout<<"Nerest Foodbank\t"<<"# of Address"<<"\t %Address"<<endl;
		cout<<"---------------\t"<<"------------"<<"\t---------"<<endl;
		cout<<"0.00 - 1.00 \t"<<countArggResult.range1<<"\t\t %"<<countArggResult.preset1<<endl;
		cout<<"1.00 - 2.00 \t"<<countArggResult.range2<<"\t\t %"<<countArggResult.preset2<<endl;
		cout<<"2.00 - 5.00 \t"<<countArggResult.range3<<"\t\t %"<<countArggResult.preset3<<endl;
		cout<<"    5.00 \t"<<countArggResult.range4<<"\t\t %"<<countArggResult.preset4<<endl;
	}


	//free virables
	delete []  countBuff;
	MPI_Type_free(&recType);
}


int main(int argc, char* argv[])
{

	if( MPI_Init(&argc, &argv) == MPI_SUCCESS )
	{
		double timeStart,timeEnd;

		int procRank,numProcs;
		MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
		timeStart = MPI_Wtime();
		MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
		MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
		processData(procRank,numProcs);
		MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
		timeEnd = MPI_Wtime();

		MPI_Finalize();
		if (procRank == 0) { /* use time on master node */
			cout<<"Runtime="<<timeEnd-timeStart<<endl;
		}


	}
	return 0;
}
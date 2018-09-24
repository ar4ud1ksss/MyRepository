/*!
	\file
	Program is helper with laba 1.1.1
	
	Сompilation options  "-lm -fno-stack-protector"
        
	//\details
	It read file with data, process them and write result in file.
	The data in the file should be stored in the format of 36 pairs of numbers of current and voltage.
	If the data file has error, then program print you line of error in the file.
*/

#include <stdio.h>
#include <math.h> 
#include <assert.h>			
					/// Number all measurements  

const int   NMeas         = 36;				
					/// Number measurments for each length                    
const int   Num           = 12;
 					/// Number (13-24) measurments for 2nd length                                
const int   NMeas2        = 12;               		
					/// Number (25-36) measurments for 3ed length 
const int   NMeas3        = 24;               		
					/// Length of wire №1 (cm) 
const int   Length1       = 20;                  		
					/// Length of wire №2 (cm) 
const int   Length2       = 30;                		
					/// Length of wire №3 (cm) 
const int   Length3       = 50;                 		
					/// Max value ampermeter (mA) 
const float   MaxI        = 300;               		
					/// Max value voltmeter (mV) 
const float   MaxU        = 600;                        
					/// Error voltmeter (mV)  
const float ErrorU        = 1.2;             		
					/// Error ampermeter (mA) 
const float ErrorI        = 1.2;             		
					/// Resistance voltmeter (Om) 
const int   RU            = 4000;                		
					/// Error square (10^(-3)cm^2)
const float ErrorD        = 0.01;                  	
					/// Error length (cm)
const float ErrorL        = 0.1;             		
					/// Number value of diameter
const int   NumD          = 10;				
					/// Number Pi
const float Pi            = 3.14;
					/// Max value resistance
const float MaxResistance = 10;
					/// Max value error
const float MaxError      = 1;
                                                                                           

int  Count ();
void Avg (float sumUI, float sumII, float* resistanceAvg);
void RandomError (float resistanceAvg, float sumUU, float sumII, float* errorRrandom);
void System (float maxU, float maxI, float resistanceAvg, float* errorRsystem);
void TotalErrorR (float errorRrandom, float errorRsystem ,float* errorR);
void TotalR (float resistanceAvg, float* rTotal, float* ro, float S, int length);
void TotalErrorRo (float ro, float errorR, float rTtotal, float* errorRo, float error, int length);
void Sort(float U[], float I[]);
void Diameter(float* s, float* errorS, FILE* inputDiametre);
void Meassage (int correct);
void Processing (float U[], float I[], float numExp, float* sumUI, float* sumUU, float* sumII);
void Read (float U[], float I[], int size, FILE* input);
int  Checkfile (float U[], float I[]);
void PrintCSV (float U[], float I[], FILE* outputCSV);
void Max (float U[], float I[], float* maxU , float* maxI, int numExp); 
void ZeroArray (float array[], int size);


	                                                                 
int main ()
	{
	printf ("#Laba 1.1.1 \n"
                "#(c) Kirill Shcherbina 2k!8\n\n");

	int correct = Count();
	assert ((correct >= -3) && (correct <= NMeas));

	Meassage (correct);
	
	}
/*!
	\brief
	Funcion that process data and print result 
	
*/
int Count()                                                                  
	{
	int size = NMeas;
	float U[size], I[size], R[size];
	ZeroArray (U, size);
	ZeroArray (I, size);
	ZeroArray (R, size);

	FILE* input          = fopen ("laba.txt",       "r");
	FILE* inputDiameter  = fopen ("diametr.txt",    "r");                                  
        FILE* output         = fopen ("Babka-labka.txt","w");                                
	FILE* outputCSV      = fopen ("grafik.csv",     "w");

	Read (U, I, size, input);
	
	if ((input == NULL) && (inputDiameter == NULL)) return (-3);
	if (input == NULL) return (-1);
	if (inputDiameter == NULL) return (-2);
	
	int check = Checkfile (U, I);

	assert (check >= 0);

	if (check > 0) return (check);

	float sumII1 = 0, sumUI1 = 0, sumUU1 = 0, maxU1 = 0, maxI1 = 0;
	float sumII2 = 0, sumUI2 = 0, sumUU2 = 0, maxU2 = 0, maxI2 = 0;
	float sumII3 = 0, sumUI3 = 0, sumUU3 = 0, maxU3 = 0, maxI3 = 0;
	int numExp1 = 0, numExp2 = NMeas2, numExp3 = NMeas3;

	Processing (U, I, numExp1, &sumUI1, &sumUU1, &sumII1);	
	Processing (U, I, numExp2, &sumUI2, &sumUU2, &sumII2);
	Processing (U, I, numExp3, &sumUI3, &sumUU3, &sumII3);
	
	Max (U, I, &maxU1, &maxI1, numExp1);
	Max (U, I, &maxU2, &maxI2, numExp2);
	Max (U, I, &maxU3, &maxI3, numExp3);	
	
	float s = 0, errorS = 0;
	Diameter (&s, &errorS, inputDiameter);	

	float resistanceAvg1 = 0, resistanceAvg2 = 0, resistanceAvg3 = 0;
	Avg (sumUI1, sumII1, &resistanceAvg1);
	Avg (sumUI2, sumII2, &resistanceAvg2);
	Avg (sumUI3, sumII3, &resistanceAvg3);
        	
	float errorRrandom1 = 0, errorRrandom2 = 0, errorRrandom3 = 0;
	RandomError (resistanceAvg1, sumUU1, sumII1, &errorRrandom1);
	RandomError (resistanceAvg2, sumUU2, sumII2, &errorRrandom2);
        RandomError (resistanceAvg3, sumUU3, sumII3, &errorRrandom3); 
	
	float errorRsystem1 = 0, errorRsystem2 = 0, errorRsystem3 = 0;	
	System (maxU1, maxI1, resistanceAvg1, &errorRsystem1);
	System (maxU2, maxI2, resistanceAvg2, &errorRsystem2);
	System (maxU3, maxI3, resistanceAvg3, &errorRsystem3);	

	float errorR1 = 0, errorR2 = 0, errorR3 = 0;
	TotalErrorR (errorRrandom1, errorRsystem1 ,&errorR1);
	TotalErrorR (errorRrandom2, errorRsystem2 ,&errorR2);
	TotalErrorR (errorRrandom3, errorRsystem3 ,&errorR3);
	
	float rTotal1 = 0, ro1 = 0, rTotal2 = 0, ro2 = 0, rTotal3 = 0, ro3 = 0;	
	TotalR (resistanceAvg1, &rTotal1, &ro1, s, Length1);
	TotalR (resistanceAvg2, &rTotal2, &ro2, s, Length2);
	TotalR (resistanceAvg3, &rTotal3, &ro3, s, Length3); 	
	
	float errorRo1 = 0, errorRo2 = 0, errorRo3 = 0;
	TotalErrorRo (ro1, errorR1, rTotal1, &errorRo1, errorS, Length1);
	TotalErrorRo (ro2, errorR2, rTotal2, &errorRo2, errorS, Length2);
	TotalErrorRo (ro3, errorR3, rTotal3, &errorRo3, errorS, Length3);

	float roTotal = (ro1 + ro2 + ro3)/3;
	float errorRoTotal = (errorRo1 + errorRo2 + errorRo3)/3;


	
	fprintf (output,"\nFor length=%dcm\t\t\tFor length=%dcm\t\t\tFor length=%dcm\n\n", Length1, Length2, Length3);
	fprintf (output,"Ravg=%5.3f Om\t\t\tRavg=%5.3f Om\t\t\tRavg=%5.3f Om\n", resistanceAvg1, resistanceAvg2, resistanceAvg3);
	fprintf (output,"Rtot=%5.3f Om\t\t\tRtot=%5.3f Om\t\t\tRtot=%5.3f Om\n", rTotal1, rTotal2, rTotal3);
	fprintf (output,"ErrorRand=%5.3f Om\t\tErrorRand=%5.3f Om\t\tErrorRand=%5.3f Om\n", errorRrandom1, errorRrandom2 , errorRrandom3);
	fprintf (output,"ErrorSys=%5.3f Om\t\tErrorSys=%5.3f Om\t\tErrorSys=%5.3f Om\n", errorRsystem1, errorRsystem2 , errorRsystem3);
	fprintf (output,"ErrorTot=%5.3f Om\t\tErrorTot=%5.3f Om\t\tErrorTot=%5.3f Om\n", errorR1, errorR2 , errorR3);
	fprintf (output,"\n\n\n\t\t\t\tl,cm     Ro,10^(-4)Om*cm     ErrorRo,10^(-4)Om*cm");
	fprintf (output,"\n\t\t\t     1.  %d         %5.2f\t\t%5.2f", Length1, ro1, errorRo1);
	fprintf (output,"\n\t\t\t     2.  %d         %5.2f\t\t%5.2f", Length2, ro2, errorRo2);
	fprintf (output,"\n\t\t\t     3.  %d         %5.2f\t\t%5.2f", Length3, ro3, errorRo3);
	fprintf (output,"\n\nTotal value resistivity (%3.2f+-%3.2f) * 10^(-4)Om * cm", roTotal, errorRoTotal); 
	
	Sort(U,I);
	
	PrintCSV (U, I, outputCSV);
	
	fclose (input);
	fclose (output);
        fclose (inputDiameter);
        fclose (outputCSV);
	
	
	return (0);                                                          
	}

//=============================================================================

void ZeroArray (float array[], int size)
	{
	for (int i = 0; i < size; i++)
        	{
        	assert (0 <= i && i < size);

        	array[i] = 0;
		}
	}

//=============================================================================

void Read (float U[], float I[], int size, FILE* input)
	{
	for (int line = 0; line < size; line++)
		{
		assert ((0 <= line) && (line < size));

		fscanf (input, "%f %f", &U[line], &I[line]);
		}
	}

//=============================================================================

int Checkfile (float U[], float I[])
	{	
	for (int line = 0; line < NMeas; line++)
		{
		assert ((0 <= line) && (line < NMeas));

		if ((U[line] <= 0) || (I[line] <= 0) || (U[line] > MaxU) || (I[line] > MaxI))
			{
			return (line + 1);
			}
        	}
		return (0);
	}

//=============================================================================

/*!
	\brief
	Funcion that count average resistance of wire 
	\details Funcion uses the lest squares method for obtain average resistance
	\param sumUI 
	ammount of multiplication current and voltage
	\param sumII 
	ammount of current^2
	\param sumUI 
	Averadge value resistance 
*/
void Avg (float sumUI, float sumII, float* resistanceAvg)
	{
	*resistanceAvg = sumUI / sumII;

	assert ((0 < *resistanceAvg) && (*resistanceAvg < MaxResistance)); 
	}

//=============================================================================

/*!
	\brief
	Funcion that count random error for resistance of wire  
	\param sumUU 
	ammount of voltage^2
	\param sumII 
	ammount of current^2
	\param Ravg 
	Averadge value resistance
	\param errorRrandom
	Value of random error for resistance 
*/
void RandomError (float resistanceAvg, float sumUU,float sumII, float* errorRrandom)
	{
	float sqrtN = sqrt(Num);
	*errorRrandom = (sqrt((sumUU / sumII) - (resistanceAvg * resistanceAvg))) / sqrtN;

	assert ((0 < *errorRrandom) && (*errorRrandom < MaxError));
	}

//=============================================================================

/*!
	\brief
	Funcion that count systematic error for resistance of wire
	\param maxU 
	Max value of voltage 
	\param maxI 
	Max value of current 
	\param Ravg 
	Averadge value resistance
	\param errorRsystem
	Value of systematic error for resistance   
*/
void System (float maxU, float maxI, float resistanceAvg, float* errorRsystem)
	{
	assert ((maxU > 0) && (maxI > 0));	
	*errorRsystem = resistanceAvg * sqrt((ErrorU / maxU) * (ErrorU / maxU) + (ErrorI / maxI) * (ErrorI / maxI));
	
	assert ((0 < *errorRsystem) && (*errorRsystem < MaxError));
	}

//=============================================================================

/*!
	\brief
	Funcion that count total error for resistance of wire
	\param errorRrandom 
	Value of random error for resistance 
	\param errorRsystem
	Value of systematic error for resistance 
	\param errorR 
	Total value error for resistance
*/
void TotalErrorR (float errorRrandom, float errorRsystem ,float* errorR)
	{
	*errorR = sqrt(errorRrandom * errorRrandom + errorRsystem * errorRsystem);
	}

//=============================================================================

/*!
	\brief
	Funcion that count total resistance and resistivity of wire (l=20cm)  
	\param Ravg1 
	Averadge value resistance
	\param Rtotal1 
	Total resistance of wire
	\param Ro1
	Value resistivity of wire
	\param S
	Value square of wire
*/
void TotalR (float resistanceAvg, float* rTotal, float* ro, float s, int length)
	{
	*rTotal = resistanceAvg + (resistanceAvg * resistanceAvg) / RU;
	*ro = *rTotal * s / length;
	}

//=============================================================================

/*!
	\brief
	Funcion that count total error for resistivity of wire (l=20cm)  
	\param Ro1
	Value resistivity of wire
	\param errorR1
	Total value error for resistance of wire
	\param Rtotal1 
	Total value resistance of wire
	\param errorRo1
	Value error for resistivity
	\param errorS
	Value error for square of wire
*/

void TotalErrorRo (float ro, float errorR, float rTotal, float* errorRo, float errorS, int length)
	{
	*errorRo = ro * sqrt((errorR / rTotal) * (errorR / rTotal) + (errorS * errorS) + (ErrorL / length) * (ErrorL / length));
	}

//=============================================================================

/*!
	\brief
	Funcion that sorts value of current and voltage for other length 
	\details Funcion return sorted array
	\param U[]  
	Array for return sorted data of voltage
	\param I[]
	Array for return sorted data of current
*/
void Sort(float U[], float I[])
	{	
	
	float saveU = 0;
	float saveI = 0;

	for (int j = 0; j < (NMeas-1); j++)
		{
		for (int i = j+1; i < NMeas2; i++)
			{
			if (U[i]<U[j])
				{
				saveU = U[j];
				U[j] = U[i];
				U[i] = saveU;	
				}
			if (I[i]<I[j])
				{
				saveI = I[j];
				I[j] = I[i];
				I[i] = saveI;	
				}
			}
		}
	
	
	saveU = 0;
	saveI = 0;
	for (int j = NMeas2; j < (NMeas3 - 1); j++)
		{
		for (int i = j+1; i < NMeas; i++)
			{
			if (U[i] < U[j])
				{
				saveU = U[j];
				U[j] = U[i];
				U[i] = saveU;	
				}
			if (I[i] < I[j])
				{
				saveI = I[j];
				I[j] = I[i];
				I[i] = saveI;	
				}
			}
		}

	saveU = 0;
	saveI = 0;
	for (int j = NMeas3; j < (NMeas-1) ; j++)
		{
		for (int i = j+1; i < NMeas; i++)
			{
			if (U[i] < U[j])
				{
				saveU = U[j];
				U[j] = U[i];
				U[i] = saveU;	
				}
			if (I[i] < I[j])
				{
				saveI = I[j];
				I[j] = I[i];
				I[i] = saveI;	
				}
			}
		}

	}

//=============================================================================

/*!
	\brief
	Funcion that count square of wire
	\details Funcion read file with value of diameter and process data
	\param S 
	Value square of wire
	\param errorS 
	Value error for square   
*/

void Diameter (float* s, float* errorS, FILE* inputDiametre)
	{
	float D[NumD];
	float sumD = 0, dAvg = 0;
	for (int i = 1; i <= NumD; i++)
		{
		fscanf (inputDiametre,"%f", &D[i]);
		sumD += D[i];
		}
	
	dAvg = sumD / NumD;
	*s = Pi * dAvg * dAvg * 100 / 4;
	*errorS = 2 * ErrorD / dAvg;
	
	
	}
//=============================================================================

void Processing (float U[], float I[], float numExp, float* sumUI, float* sumUU, float* sumII)
	{
	for (int i = numExp; i < numExp+12; i++)                                       
		{
		*sumUI += U[i] * I[i];
		*sumII += I[i] * I[i];
		*sumUU += U[i] * U[i]; 
		}
	}

//=============================================================================

void Meassage (int correct)
	{
	if (correct > 0)
		{
		printf ("ERROR in line %d\n",correct);
		}
	
	if (correct == 0)
		{
		printf ("Result saved in file 'Babka-labka.txt'\n");
		}

	if (correct == -1)
		{
		printf ("File 'laba.txt' wasn't opend\n");
		}

	if (correct == -2)
		{
		printf ("File 'diameter.txt' wasn't opend\n");
		}

	if (correct == -3)
		{
		printf ("Files 'laba.txt' and 'diameter.txt'  weren't opend\n");
		}
	}

//=============================================================================

void PrintCSV (float U[], float I[], FILE* outputCSV)
	{
	for (int i = 0; i < NMeas; i++)
		{
		assert ((0 <= i) && (i < NMeas));
		
		fprintf (outputCSV,"%4.1f %5.2f\n", U[i], I[i]);
		if ((i == 11) || (i == 23)) fprintf (outputCSV,"\n");	
		}

	}

//=============================================================================
	
void Max (float U[], float I[], float* maxU , float* maxI, int numExp) 
	{
	for (int i = numExp; i < numExp + 12; i++)
		{
		assert ((numExp <= i) && (i < numExp + 12));
		
		if (U[i] > *maxU)
			{
			*maxU = U[i]; 
			}
		
		if (I[i] > *maxI)	
			{
			*maxI = I[i]; 
			}
		}
	}




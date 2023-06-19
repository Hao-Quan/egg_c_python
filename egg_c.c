
#include <stdio.h>
#include <conio.h>
#include <string.h>
/*#include <iostream.h> */


#include<stdlib.h>


/*#include<malloc.h>*/
#include<math.h>

#include <sys/stat.h>   // stat
#include <stdbool.h>  

// 1 col image
//#define NI 600/*  unita' di input e di output*/
//#define STREAM 600

//// 4 col image
//#define NI 2000
//#define STREAM 2000

// 1 col video
//#define NI 1000
//#define STREAM 1000

// 4 col video
#define NI 2000
#define STREAM 2000



#define NM 0/* memoria */

#define NY 10/*  unita'  HIDDEN  (pesi) */
/******IL NUMERO DEI PESI PUO' ESSERE VARIATO ***********/



 /* numero totale campioni */
#define PAUSA 400/*numero di epoche prima di ogni pausa */
#define EPSILON  .1 /* tasso di apprendimento */
#define ALPHA  .00001/*tasso di dimenticanza */
#define MAXNI 5000
#define MAXNY 50
#define MAXNM 0
#define DELTA1 0.2
#define DELTA2 0.02



char buf[100], filein[100];


//int serie[STREAM*150],ep;

int i, f, h, k, j, rad, c, cs, ne, ni, nx, nm, ny, epoca, numtot, ping, nping,
countmd, countmi, pausa,
conta, g, ind[STREAM], out[STREAM][MAXNY],
cp[STREAM][MAXNY],
count[STREAM][MAXNY], s, mx, cmin,
add, bongo[STREAM], num, posr[MAXNI + MAXNM][MAXNY], neg[MAXNI + MAXNM][MAXNY],
sumw, conf[MAXNI + MAXNM], memo[STREAM], eqp, eqn,
conty, dmass, dmass1, dmass2, dmass3, cz[STREAM];

float   a[STREAM][MAXNY], t, epsilon, alpha, y[STREAM], sum[STREAM][MAXNY], meanx[STREAM],
w1[MAXNY], w0[MAXNY * 2], xtot, wtot, x[STREAM][MAXNY],
dmin[STREAM], r;



/*****************variabili z-score *************/

int eqz, zint[STREAM][MAXNY],
contaz[STREAM];
float add0, mean[STREAM], sigma[MAXNY], add1[MAXNY], add2,
z[STREAM][MAXNY], delta1, delta2;

int flg = -1;

/*********************************************************/
FILE* result; /*file di uscita dalla rete */
FILE* result1;
FILE* pesi1;

FILE* datiin; /* file di interi sorgente o output desiderato */
FILE* datiout;  /* file di reali usciti dal canale e da demodulare */
FILE* contazs;  /* file dei contatori di frequenza zscore */
FILE* zscore;    /* file  zscore */
FILE* dynsys;
FILE* pesi;
FILE* zreali;  /* file di reali da cui si trae zecore per soglia */
/************************************************************************************************/

/************************************************************************************************/
void Output()
{


	/**********RICERCA ATTRATTORI *******************/


/***********************Z-SCORE ********************************/


	unsigned int zeros = 0;
	unsigned int nzeros = 0;

	for (i = 0; i < nping; i++)

	{

		printf("\nScore reali strato %d\n\n", i);

		add0 = 0;


		for (g = 0; g < ny; g++)

		{
			//printf(" cp  %d",cp[i][g]);
			zint[i][g] = 0;
			add0 = add0 + cp[i][g];

		}

		//	gets(buf);
		mean[i] = add0 / NY;
		//	printf(" \ni %d",i);

			// printf("   %f",mean[i]);

		add2 = 0;



		for (g = 0; g < ny; g++)
		{
			add1[g] = 0;
			add1[g] = (cp[i][g] - mean[i]) * (cp[i][g] - mean[i]);



			add2 = add2 + add1[g];

			//printf(" %f",add1[g]);
		}

		sigma[i] = sqrt(add2);

		//printf(" %f",sigma[i]);



		for (g = 0; g < ny; g++)
		{


			z[i][g] = (cp[i][g] - mean[i]) / sigma[i];

			//printf("\n i %d",i);

   ////////////////////////////////////////////////////////////////////////

			printf(" %f ", z[i][g]);

			//getch();


			fprintf(zreali, "%f ", z[i][g]);

			//////////////////////////////////////////////////////////////////////////
			if (z[i][g] < -delta1)
				zint[i][g] = 0;
			else

				if (z[i][g] > delta2)

					zint[i][g] = 1;
				else
					zint[i][g] = 8;

			// printf(" %d ",zint[i][g]);



		}

		printf("\n\n\n");
		fprintf(zreali, "\n\n\n");

	}
	//printf("\ping %d",ping);
	// gets(buf);






	for (i = 0; i < nping; i++)
	{
		memo[i] = 0;

		fprintf(zreali, "\n\n");
		fprintf(zscore, "\n\n");

		printf("\n\n");

		printf("\nz-score strato %d\n\n", i);
		for (g = 0; g < ny; g++)

		{


			if (zint[i][g] == 0)
				zeros++;
			else
				nzeros++;


			fprintf(zscore, "%d ", zint[i][g]);




			printf(" %d ", zint[i][g]);

		}

		printf("\n\n");



	}





	printf("\n\n\nNumero di elementi zero:     %d", zeros);
	printf("\n\nNumero di elementi non-zero: %d", nzeros);
	printf("\n\n\n");

	//gets(buf);     

	fprintf(zscore, "\n\n\nNumero di elementi zero:     %d", zeros);
	fprintf(zscore, "\n\nNumero di elementi non-zero: %d", nzeros);









}

/************************************************************************************************/


  /*********************************************************************/
void Esegui()
{

	for (i = 0; i < nping; i++)
	{

		dmin[i] = 0;


		for (k = 0; k < ny; k++)

		{
			sum[i][k] = 0.0;

			count[i][k] = 0;
		}
	}




	/*INIZIA IL CICLO SUI VETTORI DI INPUT */

	for (i = 0; i < nping; i++)

	{
		meanx[i] = 0;
		for (h = 0; h < nx; h++)
		{

			meanx[i] = meanx[i] + x[i][h];

			for (k = 0; k < ny; k++)
			{


				a[h][k] = fabs(w1[k] - x[i][h]);

				sum[i][k] = sum[i][k] + a[h][k];


				//	printf("\n a %f",a[h][k]);
				 //	printf("\n sum %f",sum[i][k]);
				 //	  gets(buf);


			}


			meanx[i] = meanx[i] / nx; /* intanto calcolo la media sul record */

			//	printf("\nmean %f", meanx[i]);

		}


		for (k = 0; k < ny; k++)
		{

			dmin[i] = sum[i][k];
		}


		for (k = 0; k < ny; k++)
		{

			//		printf(" dmin %f",dmin[i]);
			//		printf(" sum %f",sum[i][k]);

			if (dmin[i] >= sum[i][k])
			{
				dmin[i] = sum[i][k];
			}



		}



	} /* chiude ciclo i<=nping */


	for (i = 0; i < nping; i++)
	{

		for (k = 0; k < ny; k++)
		{

			out[i][k] = 0;
		}
	}



	for (i = 0; i < nping; i++)
	{

		for (k = 0; k < ny; k++)
		{

			if (dmin[i] == sum[i][k])
			{
				//printf("g %d ",k);
				ind[i] = k;
				// printf("\n indice neurone vincente %d",k);

				out[i][k] = 1;
				/*	gets(buf); */

	  /******** SE IL NEURONE W E' VINCENTE OUT = 1 ***********/

  /************************************  NEIGHBOURHOOD  ********************************************/
	  /*VINCONO ANCHE I W DELL'INTORNO DI RAGGIO rad *****


				  if (g >= rad)
						  {out[i][g-rad][ping]=1;}
				  if (g <= ny-rad)
						  {out[i][g+rad][ping]=1;}

  /***************************************************************************************************************************/


			}


		}


		fprintf(dynsys, "%d  ", ind[i]);

		//printf("%d",ind[i]);
	}  /* chiude il ciclo su tutti gli x di input */

	printf("\n nuovo strato");
	fprintf(pesi, "\n nuovo strato  ");
	for (g = 0; g < ny; g++)
	{


		printf(" w1 %f", w1[g]);
		fprintf(pesi, " %f  ", w1[g]);

	}











	/****************************************************************/


	for (i = 0; i < nping; i++)
	{

		// printf(" \ni %d", i);
	 //printf(" \nout %d", out[i][g]);

		for (g = 0; g < ny; g++)
		{



			if (out[i][g] == 1)
				/*SE IL NEURONE W E' VINCENTE O APPARTIENE AL VICINATO MODIFICALO
				OPPORTUNAMENTE  (epsilon), ALTRIMENTI PUNISCILO (alpha) */
			{

				cp[i][g] = cp[i][g] + 1;
				//  printf(" \ni %d", i);
				 // printf(" \ncp %d", cp[i][g]);
				/* COSCIENZA :NON FA VINCERE TROPPO LO STESSO NEURONE -
				SI PUO' MODIFICARE QUANTE VOLTE FARLO VINCERE   */


				count[i][g] = count[i][g] + 1;
				if ((count[i][g] == 1))
				{

					/*??????????????????????? quale x va considerato ? */

					w1[g] = w1[g] + epsilon * (w1[g] - meanx[i]);

					//printf("\n w1 %f",w1[g]);
					//fprintf(pesi,"\n %f  ",w1[g]);



					/* EPSILON ED ALFA DEVONO DECRESCERE NEL TEMPO MA MOLTO LENTAMENTE */


					/* SE IL NEURONE VINCENTE FA PARTE DEL MIDDAMBOLO E L'OUTPUT E' QUELLO
					  GIUSTO PREMIALO, ALTRIMENTI PUNISCILO  */




				}


				/*SE VINCE TROPPO AZZERA IL CONTATORE E PUNISCILO UN PO'   */

				else
				{
					count[i][g] = 0;
					w1[g] = w1[g] - alpha * (w1[g] - meanx[i]);


				}




			}    /*chiude out = 1   */

			else /* SE IL NEURONE NON E' FRA QUELLI VINCENTI PUNISCILO  */


			{

				/* printf(" \npassa per knonvincente= %d",k);  */


				w1[g] = w1[g] - alpha * (w1[g] - meanx[i]);


			}



		}



	}



	/*chiude indici */








}





/************************************************************************************************/




float Apprendi()
{


	nping = numtot / nx;


	zscore = fopen("zeta.txt", "w+");

	/* LETTURA DEGLI INPUT PER IL TRAINING: INPUT DAL CANALE (FLOAT) */

	printf("\n\n\nNome file:     %s", filein);
	fprintf(zscore, "\n\n\nNome file:     %s", filein);


	printf("\n\n\nNumero di neuroni di input:     %d", ni);
	fprintf(zscore, "\n\n\nNumero di neuroni di input:     %d", ni);

	printf("\n\n\nNumero di campioni prelevati:     %d", numtot);
	fprintf(zscore, "\n\n\nNumero di campioni prelevati:     %d", numtot);

	printf("\n\n\nNumero di neuroni dello strato competitivo:     %d", ny);
	fprintf(zscore, "\n\n\nNumero di neuroni dello strato competitivo:     %d", ny);

	printf("\n\n\nTasso di apprendimento:     %f", epsilon);
	fprintf(zscore, "\n\n\nTasso di apprendimento:     %f", epsilon);

	printf("\n\n\nDelta1:     %f", delta1);
	fprintf(zscore, "\n\n\nDelta1:     %f", delta1);
	printf("\n\n\nDelta2:     %f", delta2);
	fprintf(zscore, "\n\n\nDelta2:     %f", delta2);


	printf("\n\n\n");
	fprintf(zscore, "\n\n\n");



	zreali = fopen("zreali.txt", "w+");


	datiin = fopen(filein, "r");

	do
	{
		fscanf(datiin, "%f", &y[i]);
		// printf("\n  %f",y[i]);


		i = i + 1;
	} while (i < numtot);

	//	gets(buf);

	fclose(datiin);


	h = 0;
	for (i = 0; i < nping; i++)

	{
		for (j = 0; j < nx; j++)

		{
			x[i][j] = y[h];
			h = h + 1;
			// divisione per 1000 solo per files .wav
				//printf("%f",x[i][j]);
				//gets(buf);

		}

	}












	epoca = 1;




	do

	{


		/* CICLO SULLE EPOCHE */





		Esegui();





		printf("\nEpoca %d\n", epoca);




		/* RAGGIO DELLA NEIGHBORHOOD - DEVE SCENDERE IL PIU' LENTAMENTE POSSIBILE */
			   /*	printf(" rad %d ",rad);
				   if(epoca%10==0)   /*per rallentare la discesa del raggio
					   {rad=rad - 1;}   */

					   /*alpha= alpha - .0001;
					   epsilon = epsilon - .0001;   */

		fprintf(dynsys, "\n");
		epoca = epoca + 1;



		if (epoca % pausa == 0)

		{

			/*printf("\nAncora? (s/n)\n");


			scanf("%c", &c);
			if (c == 'n')
			{
				epoca = -1;

				printf("\n\n\n");
			}*/

			epoca = -1;

		}



	}

	while (epoca != -1);  /**********CHIUDE CICLO EPOCHE ****************/



	Output();

	return(0);

	fclose(zscore);
	fclose(zreali);

}

/***********************************************/
 /***********************************************/
//void Setta(char *argv)
void Setta(char *argv)

{
	//printf("\nNome del file ?  ");
	/*scanf("%s", &filein);*/
	//scanf("%s", "29_beta_P7.txt");

	//filein = "29_beta_P7.txt";
	/*strncpy(filein, "29_beta_P7.txt", sizeof("29_beta_P7.txt"));*/
	strncpy(filein, argv, sizeof(argv));
	ni = NI;

	ny = NY;

	rad = NY / 2;
	nm = NM;
	alpha = ALPHA;
	epsilon = EPSILON;

	pausa = PAUSA;

	numtot = STREAM;
	delta1 = DELTA1;
	delta2 = DELTA2;


	nx = (ni + nm);




	datiin = fopen(filein, "r");
	h = 0;
	do
	{

		fscanf(datiin, "%f", &w0[k]);

		if (k % 2 == 0)
		{
			w1[h] = w0[k];
			h = h + 1;
		}


		k = k + 1;


	} while (k < ny * 2);

	//	gets(buf);
	fclose(datiin);

	/*for(h=0;h<ny;h++)
		{
		  printf("\n  %f",w0[h]);
		 printf(" new %f",w1[h]);
		}   */





		/*pesi1 = fopen("c:/baars/itsom/rita/wwdec.txt","r");
		for (k=0;k<ny;k++)
		{

				fscanf(pesi1,"%f",&w1[k]);

				 printf("\n w1 %f",w1[k]);



			}

			fclose(pesi1); */



			/* Setto i pesi */


		   /* for(k=0; k<ny; k++)
		  {
		  w1[k+1]=w1[k]+ .02;
		  if(i%2==0)
		  w1[i]=-w1[i];

		   //printf("w1 %f",w1[i]);
		  }

		   */





	Apprendi();




}


//void ff(char* line)
//{
//	strncpy(filein, line, sizeof(line));
//	ni = NI;
//
//	ny = NY;
//
//	rad = NY / 2;
//	nm = NM;
//	alpha = ALPHA;
//	epsilon = EPSILON;
//
//	pausa = PAUSA;
//
//	numtot = STREAM;
//	delta1 = DELTA1;
//	delta2 = DELTA2;
//
//
//	nx = (ni + nm);
//
//	datiin = fopen(filein, "r");
//	h = 0;
//	do
//	{
//
//		fscanf(datiin, "%f", &w0[k]);
//
//		if (k % 2 == 0)
//		{
//			w1[h] = w0[k];
//			h = h + 1;
//		}
//
//
//		k = k + 1;
//
//
//	} while (k < ny * 2);
//
//	//	gets(buf);
//	
//}

//int file_exists(char* filename) {
//	struct stat   buffer;
//	if (stat(filename, &buffer) == 0) {
//		return 1;
//	}
//	else {
//		return 0;
//	}
//}

/************************************************************************************************/
int main(int argc, char **argv)
{
	pesi = fopen("pesi.txt", "w+");
	dynsys = fopen("dyndec.txt", "w+");

	printf("%s", argv[1]);

	char line[100];
	//strncpy(line, "C:\\Users\\haoqu\\source\\repos\\egg_c\\29_beta_P7.txt", sizeof("C:\\Users\\haoqu\\source\\repos\\egg_c\\29_beta_P7.txt"));
	//strncpy(line, "C:/Users/haoqu/source/repos/egg_c/29_beta_P7.txt", sizeof("C:/Users/haoqu/source/repos/egg_c/29_beta_P7.txt"));
	strncpy(line, "29_beta_P7.txt", sizeof("29_beta_P7.txt"));

	//Setta(line);

	/************SETTA()*******************************************/

	/*strncpy(filein, line, sizeof(line));*/
	int leng = sizeof(argv[1]);

	strncpy(filein, argv[1], strlen(argv[1]));
	ni = NI;


	ny = NY;

	rad = NY / 2;
	nm = NM;
	alpha = ALPHA;
	epsilon = EPSILON;

	pausa = PAUSA;

	numtot = STREAM;
	delta1 = DELTA1;
	delta2 = DELTA2;


	nx = (ni + nm);

	datiin = fopen(filein, "r");
	h = 0;
	do
	{

		fscanf(datiin, "%f", &w0[k]);

		if (k % 2 == 0)
		{
			w1[h] = w0[k];
			h = h + 1;
		}


		k = k + 1;


	} while (k < ny * 2);
	Apprendi();


	/********************************************************/
	fprintf(dynsys, "\n\n\nNome file:     %s", filein);


	fprintf(dynsys, "\n\n\nNumero di neuroni di input:     %d", ni);


	fprintf(dynsys, "\n\n\nNumero di campioni prelevati:     %d", numtot);

	fprintf(dynsys, "\n\n\nNumero di neuroni dello strato competitivo:     %d", ny);


	fprintf(dynsys, "\n\n\nTasso di apprendimento:     %f", epsilon);

	fprintf(dynsys, "\n\n\nDelta1:     %f", delta1);
	fprintf(dynsys, "\n\n\nDelta2:     %f", delta2);



	fprintf(dynsys, "\n\n\n");

	fprintf(pesi, "\n\n\nNome file:     %s", filein);


	fprintf(pesi, "\n\n\nNumero di neuroni di input:     %d", ni);


	fprintf(pesi, "\n\n\nNumero di campioni prelevati:     %d", numtot);

	fprintf(pesi, "\n\n\nNumero di neuroni dello strato competitivo:     %d", ny);


	fprintf(pesi, "\n\n\nTasso di apprendimento:     %f", epsilon);

	fprintf(pesi, "\n\n\nDelta1:     %f", delta1);
	fprintf(pesi, "\n\n\nDelta2:     %f", delta2);


	fprintf(pesi, "\n\n\n");


	fclose(dynsys);
	fclose(pesi);



	/*dynsys=fopen("c:dyndec.txt","r");

	   for(ep=0;ep<10;ep++)
						{

		   printf("\n");
			 for(i=0;i<nping;i++)
							 {


		   fscanf(dynsys,"%d  ",&serie[i]);
		   printf("  %d",serie[i]);

							   }

						 }

		   gets(buf);     */




}   /* chiude main  */



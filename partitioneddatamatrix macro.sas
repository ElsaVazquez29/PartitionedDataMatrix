/*****************************************************************/
/**	SAS Code: %MVINTEGRATION Macro   							**/
/** Programmer: Katherine Cai									**/
/** Description: Calculates multivariate normal probabilities,  **/
/**		macro code adapted from Genz and Bretz                  **/
/** Notes: 														**/
/*****************************************************************/


%macro MVIntegration(reflib=);
/*
   SAS/IML program for the calculation of multivariate normal probabilities. The code uses the RANUNI function
   for the generation of uniform random variables. The program evaluates the multivariate normal integral by
   applying randomised lattice rule on the transfomed integral as described by Genz (1992, 1993).
   For the evaluation of singular integrals the method follow the representation of Genz and Kwong (2000).
   Further more, variable prorization and anthitec sampling is used. The program computes
   multivariate normal probabilities for positive semi-definite covariance matrices until diemnsion 100.

   Author : ALAN GENZ and FRANK BRETZ
   Contact: bretz@ifgb.uni-hannover.de
   Date   : 21.06.00 - start version

   Input  : N     : dimension of the problem (scalar)
            LOWER : lower integration limits (N-rowvector)
            UPPER : upper integration limits (N-rowvector)
            INFIN : limit flags (N-rowvector): if INFIN(I) < 0, Ith limit is (-infinity, infinity)
                                               if INFIN(I) = 0, Ith limit is (-infinity, UPPER(I)]
                                               if INFIN(I) = 1, Ith limit is [ LOWER(I), infinity)
                                               if INFIN(I) = 2, Ith limit is [ LOWER(I), UPPER(I)]
            COVAR : positive semi-definie covariance matrix (N*N-matrix)
            MAXPTS: maximum nuber of function values (scalar)
            ABSEPS: absolute error tolerance (scalar)
            RELEPS: relative error tolerance (scalar)

   Output : ERROR : estimated absolute error, with 99% confidence level
            VALUE : estimated integral value
            NEVALS: number of evaluations
            INFORM: information parameter: if INFORM = 0 then normal completion with ERROR < EPS
                                           if INFORM = 1 then completion with ERROR > EPS
                                           if INFORM = 2 then N > 100 or N < 1
                                           if INFORM = 3 then one INFIN(I) > 2 or A(I) > B(I)
                                           if INFORM = 4 then COVAR not positive semidefinite
*/


/*

Example call:
*************

The statements at the end of the program


  N = 5;

  LOWER = J(1,N,-2);
  UPPER = J(1,N,2);
  INFIN = J(1,N,2);

  COVAR = {   1    0 0    0 0,
              0    1 0    0 0,
              0    0 1    0 0,
              0    0 0    1 0,
              0    0 0    0 1};

  MAXPTS = 300;*2000*N*N*N;
  ABSEPS = .0001;
  RELEPS = 0;

  RUN MVN_DIST( N, LOWER, UPPER, INFIN, COVAR, MAXPTS, ABSEPS, RELEPS,  ERROR, VALUE, NEVALS, INFORM );

  PRINT ERROR VALUE NEVALS INFORM;


lead to the following output:

    ERROR     VALUE    NEVALS    INFORM
0.0000756 0.9030463     27040         0

*/



/*OPTIONS STIMER NOCENTER NOSOURCE;*/

PROC IML SYMSIZE = 250000;

NL = 100;

START MVN_DIST( N, LOWER, UPPER, INFIN, COVAR, MAXPTS, ABSEPS, RELEPS,   ERROR, VALUE, NEVALS, INFORM );
      NEVALS = 0;
      RUN MVNDNT( N, COVAR, LOWER, UPPER, INFIN,   INFIS, VALUE, ERROR, INFORM );
      IF ( INFORM = 0 ) THEN DO;
         IF ( N-INFIS > 2 ) THEN RUN DKBVRC( N-INFIS-1, 0, MAXPTS, ABSEPS, RELEPS, ERROR, VALUE, NEVALS, INFORM );
      END;
FINISH MVN_DIST;

START MVN_DFN( N, W ) GLOBAL( COVARS, DONE, EONE, INFI, A, B, Y );
      VALUE = 1;
      INFA = 0;
      INFB = 0;
      IK = 1;
      DO I = 1 TO N+1;
         VSUM = 0;
         IF ( IK > 1 ) THEN VSUM = VSUM + SUM( COVARS[I,1:IK-1]#T(Y[1:IK-1]));
         IF ( INFI[I] ^= 0 ) THEN DO;
            IF ( INFA = 1 ) THEN AI = MAX( AI, A[I] - VSUM );
            ELSE DO;
               AI = A[I] - VSUM;
               INFA = 1;
            END;
         END;
         IF ( INFI[I] ^= 1 ) THEN DO;
            IF ( INFB = 1 ) THEN BI = MIN( BI, B[I] - VSUM );
            ELSE DO;
               BI = B[I] - VSUM;
               INFB = 1;
            END;
         END;
         IF ( I < N + 1 & IK < N + 1 ) THEN AAA = COVARS[I+1,IK+1];
         ELSE AAA = 0;
         IF ( I = N+1 | AAA > 0 ) THEN DO;
            IF ( I = 1 ) THEN DO;
               DI = DONE;
               EI = EONE;
            END;
            ELSE DO;
               DI = 0;
               EI = 1;
               J = 2*INFA+INFB-1;
               IF ( J >= 0 ) THEN DO;
                  IF ( J ^= 0 ) THEN DI = PROBNORM(AI);
                  IF ( J ^= 1 ) THEN EI = PROBNORM(BI);
               END;
               EI = MAX( EI, DI );
            END;
            IF ( DI >= EI ) THEN DO;
               VALUE = 0;
               I = N+1;
            END;
            ELSE DO;
               VALUE = VALUE*( EI - DI );
               IF ( I <= N ) THEN Y[IK] = PROBIT( DI + W[IK]*( EI - DI ) );
               IK = IK + 1;
               INFA = 0;
               INFB = 0;
            END;
         END;
      END;
      MVNDFN = VALUE;
      RETURN( MVNDFN );
FINISH MVN_DFN;


START MVNDNT( N, COVAR, LOWER, UPPER, INFIN,   INFIS, VALUE, ERROR, INFORM ) GLOBAL( COVARS, DONE, EONE, INFI, A, B, NL );
      INFORM = 0;
      IF ( N > NL | N < 1 ) THEN INFORM = 2;
      ELSE DO  I = 1 TO N;
         IF ( INFIN[I] > 2 ) THEN INFORM = 3;
         ELSE IF ( INFIN[I] = 2 & LOWER[I] > UPPER[I] ) THEN INFORM = 3;
      END;
      IF ( INFORM = 0 ) THEN RUN COVSRT( N, LOWER, UPPER, COVAR, INFIN, INFIS, INFORM );
      IF ( INFORM = 0 ) THEN DO;
         IF ( N - INFIS = 0 ) THEN DO;
            VALUE = 1;
            ERROR = 0;
         END;
         ELSE DO;
            IF ( N - INFIS = 1 ) THEN DO;
               VALUE = EONE - DONE;
               ERROR = 2E-15;
            END;
            ELSE DO;
               IF ( N - INFIS = 2 ) THEN DO;
                  IF ( ABS( COVARS[2,2] ) > 0 ) THEN DO;
                     D = SQRT( 1 + COVARS[2,1]**2 );
                     IF ( INFI[2] ^= 0 ) THEN A[2] = A[2]/D;
                     IF ( INFI[2] ^= 1 ) THEN B[2] = B[2]/D;
                     VALUE = PROBBVN( A, B, INFI, COVARS[2,1]/D );
                  END;
                  ELSE DO;
                     IF ( INFI[1] ^= 0 ) THEN DO;
                        IF ( INFI[2] ^= 0 ) THEN A[1] = MAX( A[1], A[2] );
                     END;
                     ELSE DO;
                        IF ( INFI[2] ^= 0 ) THEN A[1] = A[2];
                     END;

                     IF ( INFI[1] ^= 1 ) THEN DO;
                        IF ( INFI[2] ^= 1 ) THEN B[1] = MIN( B[1], B[2] );
                     END;
                     ELSE DO;
                        IF ( INFI[2] ^= 1 ) THEN B[1] = B[2];
                     END;

                     IF ( INFI[1] ^= INFI[2] ) THEN INFI[1] = 2;
                     RUN MVNLMS( A[1], B[1], INFI[1],   D, E );
                     VALUE = E - D;
                  END;
                  ERROR = 2E-15;
               END;
               ELSE DO;
                  VALUE = 0;
                  ERROR = 1;
               END;
            END;
         END;
      END;
      ELSE DO;
         VALUE = 0;
         ERROR = 1;
      END;
FINISH MVNDNT;

START MVNLMS( A, B, INFIN,   LOWER, UPPER );
      LOWER = 0;
      UPPER = 1;
      IF ( INFIN >= 0 ) THEN DO;
         IF ( INFIN ^= 0 ) THEN LOWER = PROBNORM(A);
         IF ( INFIN ^= 1 ) THEN UPPER = PROBNORM(B);
      END;
      UPPER = MAX( UPPER, LOWER );
FINISH MVNLMS;

START COVSRT( N, LOWER, UPPER, COVAR, INFIN,   INFIS, INFORM )
      GLOBAL( EPS, SQTWPI, COVARS, DONE, EONE, INFI, A, B, Y );

      Y = J( 1, N, . );
      INFI = INFIN;
      A = J( 1, N, 0 );
      B = J( 1, N, 0 );
      COVARS = COVAR;
      INFIS = N - SUM( SIGN( SIGN( INFI ) + 1 ) );
      DO I = 1 TO N;
         IF ( INFI[I] >= 0 ) THEN DO;
            IF ( INFI[I] ^= 0 ) THEN A[I] = LOWER[I];
            IF ( INFI[I] ^= 1 ) THEN B[I] = UPPER[I];
         END;
      END;

      IF ( INFIS < N ) THEN DO;
         DO I = N TO N-INFIS+1 BY -1;
            IF ( INFI[I] >= 0 ) THEN DO;
               DO J = 1 TO I-1;
                  IF ( INFI[J] < 0 ) THEN DO;
                     RUN RCSWP( J, I, A, B, INFI, N, COVARS );
                     J = I-1;
                  END;
               END;
            END;
         END;

         DO I = 1 TO N-INFIS;
            DEMIN = 1;
            JMIN = I;
            CVDIAG = 0;
            EPSI = I*I*EPS;
            DO J = I TO N-INFIS;
               IF ( COVARS[J,J] > EPSI ) THEN DO;
                  SUMSQ = SQRT( COVARS[J,J] );
                  VSUM = 0;
                  IF ( I > 1 ) THEN VSUM = SUM( COVARS[J,1:I-1] # T(Y[1:I-1]) );
                  AJ = ( A[J] - VSUM )/SUMSQ;
                  BJ = ( B[J] - VSUM )/SUMSQ;
                  RUN MVNLMS( AJ, BJ, INFI[J],   DD, EE );
                  IF ( DEMIN >= EE - DD ) THEN DO;
                     JMIN = J;
                     AMIN = AJ;
                     BMIN = BJ;
                     DEMIN = EE - DD;
                     CVDIAG = SUMSQ;
                  END;
               END;
            END;

            IF ( JMIN > I ) THEN RUN RCSWP( I, JMIN, A, B, INFI, N, COVARS );

            IF ( CVDIAG > 0 ) THEN DO;
               COVARS[I,I] = CVDIAG;
               DO L = I+1 TO N-INFIS;
                  COVARS[L,I] = COVARS[L,I]/CVDIAG;
                  COVARS[L,I+1:L] = COVARS[L,I+1:L] - COVARS[L,I] # T(COVARS[I+1:L,I]);
               END;

               IF ( DEMIN > EPSI ) THEN DO;
                  YL = 0;
                  YU = 0;
                  IF ( INFI[I] ^= 0 ) THEN YL = -EXP( -AMIN**2/2 )/SQTWPI;
                  IF ( INFI[I] ^= 1 ) THEN YU = -EXP( -BMIN**2/2 )/SQTWPI;
                  Y[I] = ( YU - YL )/DEMIN;
               END;
               ELSE DO;
                  IF ( INFI[I] = 0 ) THEN Y[I] = BMIN;
                  IF ( INFI[I] = 1 ) THEN Y[I] = AMIN;
                  IF ( INFI[I] = 2 ) THEN Y[I] = ( AMIN + BMIN )/2;
               END;

               COVARS[I,1:I] = COVARS[I,1:I]/CVDIAG;
               A[I] = A[I]/CVDIAG;
               B[I] = B[I]/CVDIAG;
            END;
            ELSE DO;
               IF ( COVARS[I,I] > -EPSI ) THEN DO;
                  COVARS[I:N-INFIS,I] = 0;

                  AAA = 0;
                  DO J = I-1 TO 1 BY -1;
                     IF ( ABS( COVARS[I,J] ) > EPSI ) THEN DO;
                        A[I] = A[I]/COVARS[I,J];
                        B[I] = B[I]/COVARS[I,J];
                        IF ( COVARS[I,J] < 0 ) THEN DO;
                           AA = A[I];
                           A[I] = B[I];
                           B[I] = AA;
                           IF ( INFI[I] ^= 2 ) THEN INFI[I] = 1 - INFI[I];
                        END;
                        COVARS[I,1:J] = COVARS[I,1:J]/COVARS[I,J];
                        DO L = J+1 TO I-1;
                           IF( COVARS[L,J+1] > 0 ) THEN DO;
                              DO K = I-1 TO L BY -1;

                                 AA = COVARS[K,1:K];
                                 COVARS[K,1:K] = COVARS[K+1,1:K];
                                 COVARS[K+1,1:K] = AA;

                                 AA = A[K];
                                 A[K] = A[K+1];
                                 A[K+1] = AA;

                                 AA = B[K];
                                 B[K] = B[K+1];
                                 B[K+1] = AA;

                                 M = INFI[K];
                                 INFI[K] = INFI[K+1];
                                 INFI[K+1] = M;
                              END;
                              L = I-1;
                           END;
                        END;
                        J = 1;
                        AAA = 1;
                     END;
                     IF AAA = 1 THEN;
                     ELSE COVARS[I,J] = 0;
                  END;
                  Y[I] = 0;
               END;
               ELSE DO;
                 INFORM = 4;
                 I = N-INFIS;
               END;
            END;
         END;
         IF (INFORM = 0 ) THEN RUN MVNLMS( A[1], B[1], INFI[1], DONE, EONE );
      END;
FINISH COVSRT;


START RCSWP( P, Q, A, B, INFIN, N, C );
      AA = A[P];
      A[P] = A[Q];
      A[Q] = AA;

      AA = B[P];
      B[P] = B[Q];
      B[Q] = AA;

      I = INFIN[P];
      INFIN[P] = INFIN[Q];
      INFIN[Q] = I;

      AA = C[P,P];
      C[P,P] = C[Q,Q];
      C[Q,Q] = AA;

      IF (P>1) THEN DO;
         AA = C[Q,1:P-1];
         C[Q,1:P-1] = C[P,1:P-1];
         C[P,1:P-1] = AA;
      END;

      DO I = P+1 TO Q-1;
         AA = C[I,P];
         C[I,P] = C[Q,I];
         C[Q,I] = AA;
      END;

      IF (Q<N) THEN DO;
         AA = C[Q+1:N,P];
         C[Q+1:N,P] = C[Q+1:N,Q];
         C[Q+1:N,Q] = AA;
      END;
FINISH RCSWP;

START DKBVRC( NDIM, MINVLS, MAXVLS, ABSEPS, RELEPS,   ABSERR, FINEST, INTVLS, INFORM )
      GLOBAL( PLIM, KLIM, P, C, MINSMP );

      VK     = J( 1, KLIM, . );
      INFORM = 1;
      INTVLS = 0;
      KLIMI  = KLIM;
      IF ( MINVLS >= 0 ) THEN DO;
         FINEST = 0;
         VAREST = 0;
         SAMPLS = MINSMP;
         DO I = 1 TO PLIM;
            NP = I;
            IF ( MINVLS < 2*SAMPLS*P[I] ) THEN I = PLIM;
         END;
         IF ( MINVLS >= 2*SAMPLS*P[PLIM] ) THEN SAMPLS = MINVLS/( 2*P[PLIM] );
      END;
      VALUE = J( 1, SAMPLS, . );
      EXIT  = 0;
      DO UNTIL( EXIT = 1);
         VK[1] = 1/P[NP];
         DO I = 2 TO MIN( NDIM, KLIM );
            VK[I] = MOD( C[NP, MIN(NDIM-1,KLIM-1)]*VK[I-1], 1 );
         END;
         FINVAL = 0;
         VARSQR = 0;
         DO I = 1 TO SAMPLS;
            VALUE[I] = DKSMRC( NDIM, KLIMI, P[NP], VK );
         END;
         FINVAL = VALUE[:];
         VARSQR = (VALUE[##] - VALUE[+]##2/SAMPLS) / (SAMPLS # (SAMPLS-1));
         INTVLS = INTVLS + 2*SAMPLS*P[NP];
         VARPRD = VAREST*VARSQR;
         FINEST = FINEST + ( FINVAL - FINEST )/( 1 + VARPRD );
         IF ( VARSQR > 0 ) THEN VAREST = ( 1 + VARPRD )/VARSQR;
         ABSERR = 3*SQRT( VARSQR/( 1 + VARPRD ) );
         IF ( ABSERR > MAX( ABSEPS, ABS(FINEST)*RELEPS ) ) THEN DO;
            IF ( NP < PLIM ) THEN NP = NP + 1;
            ELSE DO;
               SAMPLS = MIN( 3*SAMPLS/2, ( MAXVLS - INTVLS )/( 2*P[NP] ) );
               SAMPLS = MAX( MINSMP, SAMPLS );
            END;
            IF ( INTVLS + 2*SAMPLS*P[NP] > MAXVLS ) THEN EXIT = 1;
         END;
         ELSE DO;
            INFORM = 0;
            EXIT = 1;
         END;
      END;
FINISH DKBVRC;


START DKSMRC( NDIM, KLIM, PRIME, VK );
      X = J( 1, NDIM, . );
      SUMKRO = 0;
      NK = MIN( NDIM, KLIM );
      DO J = 1 TO NK-1;
         JP = J + RANUNI(141071)*( NK + 1 - J );
         XT = VK[J];
         VK[J] = VK[INT(JP)];
         VK[INT(JP)] = XT;
      END;
      XP = RANUNI( J(1, NDIM, 141071) );
      DIFF = NDIM - KLIM;
      DO K = 1 TO PRIME;
         IF ( DIFF>0 ) THEN DO;
            RUN DKRCHT( DIFF, X );
            DO J = 1 TO NDIM-KLIM;
               X[NK+J] = X[J];
            END;
         END;
         X[1:NK] = MOD( K*VK[1:NK], 1 );
         DO J = 1 TO NDIM;
            XT = X[J] + XP[J];
            IF ( XT > 1 ) THEN  XT = XT - 1;
            X[J] = ABS( 2*XT - 1 );
         END;
         MVNDFN = MVN_DFN( NDIM, X );
         SUMKRO = SUMKRO + ( MVNDFN - SUMKRO )/( 2*K - 1 );
         X = 1 - X;
         MVNDFN = MVN_DFN( NDIM, X );
         SUMKRO = SUMKRO + ( MVNDFN - SUMKRO )/( 2*K );
      END;
      RETURN (SUMKRO);
FINISH DKSMRC;

START DKRCHT( S, QUASI ) GLOBAL( NN, PSQT, HISUM, OLDS, MXDIM, MXHSUM, BB );

      IF ( S ^= OLDS | S < 1 ) THEN DO;
         OLDS = S;
         NN[1] = 0;
         HISUM = 0;
      END;

      I = 0;
      CRIT = 0;
      DO UNTIL( CRIT = 1 | I = HISUM + 1 );
         NN[I + 1] = NN[I + 1] + 1;
         IF ( NN[I + 1] < BB ) THEN DO;
           CRIT = 1;
           I = I - 1;
         END;
         ELSE NN[I + 1] = 0;
         I = I + 1;
      END;

      IF ( I > HISUM ) THEN DO;
         HISUM = HISUM + 1;
         IF ( HISUM > MXHSUM ) THEN HISUM = 0;
         NN[HISUM + 1] = 1;
      END;

      RN = 0;
      DO I = HISUM TO 0 BY -1;
         RN = NN[I + 1] + BB * RN;
      END;
      QUASI[1:S] = MOD( RN # PSQT[1:S], 1 );
FINISH DKRCHT;

START PROBBVN( LOWER, UPPER, INFIN, CORREL );
      IF ( INFIN[1] = 2  & INFIN[2] = 2 ) THEN
         BVN =  PROBBNRM ( LOWER[1], LOWER[2], CORREL ) - PROBBNRM ( UPPER[1], LOWER[2], CORREL )
              - PROBBNRM ( LOWER[1], UPPER[2], CORREL ) + PROBBNRM ( UPPER[1], UPPER[2], CORREL );
      ELSE IF ( INFIN[1] = 2  & INFIN[2] = 1 ) THEN
         BVN =  PROBBNRM ( -LOWER[1], -LOWER[2], CORREL ) - PROBBNRM ( -UPPER[1], -LOWER[2], CORREL );
      ELSE IF ( INFIN[1] = 1  & INFIN[2] = 2 ) THEN
         BVN =  PROBBNRM ( -LOWER[1], -LOWER[2], CORREL ) - PROBBNRM ( -LOWER[1], -UPPER[2], CORREL );
      ELSE IF ( INFIN[1] = 2  & INFIN[2] = 0 ) THEN
         BVN =  PROBBNRM ( UPPER[1], UPPER[2], CORREL ) - PROBBNRM ( LOWER[1], UPPER[2], CORREL );
      ELSE IF ( INFIN[1] = 0  & INFIN[2] = 2 ) THEN
         BVN =  PROBBNRM ( UPPER[1], UPPER[2], CORREL ) - PROBBNRM ( UPPER[1], LOWER[2], CORREL );
      ELSE IF ( INFIN[1] = 1  & INFIN[2] = 0 ) THEN BVN = PROBBNRM ( -LOWER[1],  UPPER[2], -CORREL );
      ELSE IF ( INFIN[1] = 0  & INFIN[2] = 1 ) THEN BVN = PROBBNRM (  UPPER[1], -LOWER[2], -CORREL );
      ELSE IF ( INFIN[1] = 1  & INFIN[2] = 1 ) THEN BVN = PROBBNRM ( -LOWER[1], -LOWER[2],  CORREL );
      ELSE IF ( INFIN[1] = 0  & INFIN[2] = 0 ) THEN BVN = PROBBNRM (  UPPER[1],  UPPER[2],  CORREL );
      RETURN ( BVN );
FINISH PROBBVN;



HISUM = .;
OLDS  = 0;
MXDIM = 80;
MXHSUM = 50;
BB = 2;
PSQT={1.414213562373 1.732050807569 2.236067977500 2.645751311065 3.316624790355 3.605551275464
      4.123105625618 4.358898943541 4.795831523313 5.385164807135 5.567764362830 6.082762530298
      6.403124237433 6.557438524302 6.855654600401 7.280109889281 7.681145747869 7.810249675907
      8.185352771872 8.426149773176 8.544003745318 8.888194417316 9.110433579144 9.433981132057
      9.848857801796 10.04987562112 10.14889156509 10.34408043279 10.44030650891 10.63014581273
      11.26942766958 11.44552314226 11.70469991072 11.78982612255 12.20655561573 12.28820572744
      12.52996408614 12.76714533480 12.92284798332 13.15294643797 13.37908816026 13.45362404707
      13.82027496109 13.89244398945 14.03566884762 14.10673597967 14.52583904633 14.93318452307
      15.06651917332 15.13274595042 15.26433752247 15.45962483374 15.52417469626 15.84297951775
      16.03121954188 16.21727474023 16.40121946686 16.46207763315 16.64331697709 16.76305461424
      16.82260384126 17.11724276862 17.52141546794 17.63519208855 17.69180601295 17.80449381476
      18.19340539866 18.35755975069 18.62793601020 18.68154169227 18.78829422806 18.94729532150
      19.15724406067 19.31320791583 19.46792233393 19.57038579078 19.72308292332 19.92485884517
      20.02498439450 20.22374841616};


EPS = 1E-10;
SQTWPI = 2.506628274631000502415765284811045253;

PLIM = 25;
KLIM = 20;
MINSMP = 8;
P = { 31 47 73 113 173 263 397 593 907 1361 2053 3079 4621 6947 10427 15641
      23473 35221 52837 79259 118891 178349 267523 401287 601942};
C = { 12 9 9 13 12 12 12 12 12 12 12 12 3 3 3 12 7 7 12,
      13 11 17 10 15 15 15 15 15 15 22 15 15 6 6 6 15 15 9 ,
      27 28 10 11 11 20 11 11 28 13 13 28 13 13 13 14 14 14 14 ,
      35 27 27 36 22 29 29 20 45 5 5 5 21 21 21 21 21 21 21 ,
      64 66 28 28 44 44 55 67 10 10 10 10 10 10 38 38 10 10 10 ,
      111 42 54 118 20 31 31 72 17 94 14 14 11 14 14 14 94 10 10 ,
      163 154  83 43 82 92 150 59 76 76 47 11 11 100 131 116 116 116 116 ,
      246 189 242 102 250 250 102 250 280 118 196 118 191 215 121 121 49 49 49 ,
      347 402 322 418 215 220 339 339 339 337 218 315 315 315 315 167 167 167 167 ,
      505 220 601 644 612 160 206 206 206 422 134 518 134 134 518 652 382 206 158 ,
      794 325 960 528 247 247 338 366 847 753 753 236 334 334 461 711 652 381 381 ,
      1189 888 259 1082 725 811 636 965 497 497 1490 1490 392 1291 508 508 1291 1291 508 ,
      1763 1018 1500 432 1332 2203 126 2240 1719 1284 878 1983 266 266 266 266 747 747 127 ,
      2872 3233 1534 2941 2910 393 1796 919 446 919 919 1117 103 103 103 103 103 103 103 ,
      4309 3758 4034 1963 730 642 1502 2246 3834 1511 1102 1102 1522 1522 3427 3427 3928 915 915 ,
      6610 6977 1686 3819 2314 5647 3953 3614 5115 423 423 5408 7426 423 423 487 6227 2660 6227 ,
      9861 3647 4073 2535 3430 9865 2830 9328 4320 5913 10365 8272 3706 6186 7806 7806 7806 8610 2563 ,
      10327 7582 7124 8214 9600 10271 10193 10800 9086 2365 4409 13812 5661 9344 9344 10362 9344 9344 8585 ,
      19540 19926 11582 11113 24585 8726 17218 419 4918 4918 4918 15701 17710 4037 4037 15808 11401 19398 25950 ,
      34566 9579 12654 26856 37873 38806 29501 17271 3663 10763 18955 1298 26560 17132 17132 4753 4753 8713 18624 ,
      31929  49367 10982 3527 27066 13226 56010 18911 40574 20767 20767 9686 47603 47603 11736 11736 41601 12888 32948 ,
      40701  69087 77576 64590 39397 33179 10858 38935 43129 35468 35468 2196 61518 61518 27945 70975 70975 86478 86478 ,
      103650 125480 59978 46875 77172 83021 126904 14541 56299 43636 11655 52680 88549 29804 101894 113675 48040 113675 34987 ,
      165843 90647 59925 189541 67647 74795 68365 167485 143918 74912 167289 75517 8148 172106 126159 35867 35867 35867 121694 ,
      130365 236711 110235 125699 56483 93735 234469 60549 1291 93937 245291 196061 258647 162489 176631 204895 73353 172319 28881};



libname reflib &reflib.;

RESET STORAGE=reflib.MVIntegration;
SHOW STORAGE;

STORE;

SHOW STORAGE;


QUIT;

%mend MVIntegration;







/****************************************************************************/
/**	SAS Code: %partitionedDataMatrix							           **/
/** Programmer: Elsa Vazquez Arreola									   **/
/** Description: Creates partitioned Data matrix for simultaneous modeling **/
/*****************************************************************/
%macro partitionedDataMatrix(ds=., file=, timeVar=, outVar=, predVarTD=, idVar=, alpha=0.05, predVarTI=., distr=bin, optim=NLPCG, MC=LWY);

*Check if there is a library;
%if &ds. ^=. %then %do;
	LIBNAME DS &ds.;  
	DATA mydata; 
		SET DS.&file.; 
	RUN;
%end;
%else %do;
	DATA mydata;
		SET &file.;
	RUN;
%end;

TITLE "Partitioned GMM";
PROC SORT DATA=mydata OUT=mydatasorted; 
BY &timeVar.; RUN;

*Check if there is a time independent variable;
%if &predVarTI. ^=.  %then %do;
	%let predVar = &predVarTD &predVarTI;
%end;
%else %do;
	%let predVar = &predVarTD;
%end;

*Use either the LWY or LS moment check;
%if &MC=LWY %then %do; *LWY approach;

*Obtain residuals;
%if &distr.=normal %then %do;
	proc reg data=mydatasorted NOPRINT;
	BY &timeVar.;
	MODEL &outVar. = &predVar.;
	OUTPUT OUT=outpool3 PREDICTED=mu;
	RUN;

	DATA outpool3; SET outpool3;
	wt = 1/mu; 
	rsdraw = &outVar.-mu; 
%end;
%else %if &distr.=bin %then %do;
	PROC logistic DATA=mydatasorted NOPRINT; 
	BY &timeVar.;
	MODEL &outVar. (event='1') = &predVar. / aggregate scale=none;
	OUTPUT OUT=outpool3 P=mu;
	RUN;

	DATA outpool3; SET outpool3;
	wt = mu*(1-mu);
	rsdraw = &outVar.-mu; 
%end;
%else %do;
	%put ERROR must be normal or binomial distributions; 
	%return;
%end;

PROC SORT DATA=outpool3 OUT=outpool3 ;
  BY &idVar. &timeVar.; RUN;
quit;

PROC IML;
use outpool3;                                                           
%if &predVarTI. ^=. %then %do;
	read all VARIABLES {&predVarTD. &predVarTI.} into Zmat; 
	read all var {&predVarTI.} into timeInd;
	read all var {&predVarTD.} into timeDep;
%end;
%else %do;
	read all VARIABLES {&predVarTD.} into Zmat;
%end;
read all var {wt} into wt;
read all var {rsdraw} into rsd;
read all var {&idVar.} INTO ID;
read all var {&timeVar.} INTO time;
close outpool3;


use mydata;
read all var{&predVarTD.} into timeDep2;
read all var {&predVarTI.} into timeInd2;
close mydata;



N=ncol(unique(ID)); *Number of subjects;
T = max(time); *Number of time points;
Np = ncol(Zmat); *Number of predictors, not including intercept;

NpTD2 = ncol(timeDep2);
NpTI2=ncol(timeInd2);

%if &predVarTI. ^=. %then %do;
	NpTI = ncol(timeInd);*The last NPTI columns of Zmat will be time independent;
	NpTD = ncol(timeDep); *The first NPTD columns of Zmat will be time dependent;
%end;

*Find the correlation between X (a) and Y (rsd) across the T time points;
start rho(a,rsd) global(N,T);
abm = j(N,2*T,.);
abm[,1:T] = shape(rsd,N);	* N x T - First T columns are the values of rsd sorted into T columns (by time);
abm[,T+1:2*T] = shape(a,N); * Remaining T columns are the values of a sorted into T columns;
corr = corr(abm);  
rho = corr[1:T,T+1:2*T];    * T x T - Only take the correlations between the X and Y in the T time points;
return(rho);
finish rho;


*Standard deviation for each correlation;
start stddev(a,rsd) global(N,T);
bm = shape(rsd,N);    		 * N x T;
bdev = bm-j(N,1,1)*bm[:,];   * bdev N x T,   bm[:,] Col Mean is a row vector   1 x T - each row of bm minus column means;
bdev2 = bdev#bdev;      
am = shape(a,N);   
adev = am-j(N,1,1)*am[:,];  *N x T;
adev2 = adev#adev;      
stddev = sqrt( (1/N)*t(bdev2)*adev2 );   * T x T;
return(stddev);
finish stddev;

* corrected standardization;
start stdzn(x) global(N,T);
xrows = shape(x,N);   *N x T - by shape default columns are T=nrow(x)/N;
y = xrows - xrows[:,];  *N x T - Each value minus the column mean (mean for that time point);
vcv = (1/(N-1))*t(y)*y; * T x T;
v = diag(vcv); * T x T diagonal elements of vcv;
sinv = sqrt(inv(v));
x2 = y*sinv;   *N x T;
x2  = shape(x2,N*T,1); *N*T x 1 vector of standardized values;
return(x2);
finish stdzn;

pvec = j(Np*T*T,1,.);sevec = j(Np*T*T,1,.);  * pvec   (T*T) x Np;
r4out = j(T,T*Np,.); se4out = j(T,T*Np,.); z4out = j(T,T*Np,.); p4out = j(T,T*Np,.); 

y = rsd;
y_std = stdzn(y);




DO i=1 TO Np;
x = wt#Zmat[,i]; 			 * (N*T) x 1;
x_std = stdzn(x);

*Find p-values for the correlation of X_i and Y;
r = rho(x_std,y_std);		 * T x T;
se = stddev(x_std,y_std);	 * T x T;
z = sqrt(N)*(r/se);
p = 2*(1-cdf("normal",abs(z)));  * T x T;

*Fill the corresponding T columns of each matrix with the calculated values;
r4out[,T*(i-1)+1:T*i] = r;
se4out[,T*(i-1)+1:T*i] = se;
z4out[,T*(i-1)+1:T*i] = z;
p4out[,T*(i-1)+1:T*i] = p;

DO j = 1 TO T;
p[j,j] = 1;  * Not going to test diagonal elements, set to 1;
END;
*Takes the values in p4out and se4out, across row, then down column and creates vectors;
pvec[T*T*(i-1)+1:T*T*i,1] = shape(p,T*T,1);    *(T*T) x 1;
sevec[T*T*(i-1)+1:T*T*i,1] = shape(se,T*T,1);   *(T*T) x 1;
END;


TypeVec2 = (pvec >= &alpha.*j(Np*T*T,1,1) );

Type2 = shape(TypeVec2, Np, T*T);


*Individual test approach;
x = wt; 			 * (N*T) x 1;
x_std = stdzn(x);
r = rho(x_std,y_std);		 * T x T;
se = stddev(x_std,y_std);	 * T x T;
z = sqrt(N)*(r/se);
p = 2*(1-cdf("normal",abs(z)));  * T x T;


T_wt = p>&alpha.;
T_wt = shape(T_wt,1);


TypeMtx3 = j(Np+1,T*T,.);
TypeMtx3[1,] = T_wt;
TypeMtx3[2:Np+1,] = Type2;
Type2[1,] = shape(I(T),1);


*Adjust the last NpTI rows of TypeMtx3 to account for time independent covariates;
%if &predVarTI. ^=. %then %do;
	TypeMtx3= TypeMtx3[1:(NPTD+1),]; *Subset TypeMtx3 to include only time dependent;
	TypeMtxIND = repeat(shape(I(T),1),NPTI); *Create rows of typeMtx for the time independent;
	Np = NpTD;
%end;

print typemtx3; /*typemtx3 contains all valid moments for each of the time-dependent covariates*/
/*TypeMtx3 is not affected by the number of time-independent covariates in the model*/


/* define helper functions ROW and COL */
start row(x);  /* return matrix m such that m[i,j] = i */
   return( repeat( T(1:x), 1, x ));
finish;
start col(x);  /* return matrix m such that m[i,j] = j */
   return( repeat(1:x, x) );
finish;

*Helper matrices;
rt = row(T);
ct = col(T);

PRINT RT;
PRINT CT;


*Indices of upper diagonal for removal;
upperIdx = loc(ct>rt);
lowerIdx= loc(ct<rt);


*Remove backwards conditions;
typeMtxNB = TypeMtx3;
typeMtxNB[,upperIdx] = 0;



/**************************************************************************/
/**************************************************************************/
/**************************************************************************/
/*****************  MATRIX FOR MODIFIED DATASET ***************************/
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/


MatrixCreateNew=j(T,NPTD2,1);
/*creating matrix from first time dependent covariate*/


/*creating matrix for all other time dependent covariate*/

do i=1 to T-1;
	MatrixC1=j(T,NPTD2,0);
	do j=1 to NpTD2;
      MatrixC2=shape(typeMtxNB[j+1,], T, T);
	  MatrixR=j(T-i, T-i,.);
      MatrixR[1:T-i, 1:T-i]=MatrixC2[i+1:T, 1:T-i];
      VecC3=vecdiag(MatrixR);
      MatrixC1[i+1:T, j]=vecC3;
	end;
MatrixCreateNew=MatrixCreateNew||MatrixC1;
end;

cov_id = MatrixCreateNew[+,];

keepcovar_id = loc(cov_id>0);


ncolumns=ncol(keepcovar_id);
matrixCreateNew2=matrixCreateNew[,keepcovar_id];



MatrixCreateNew3 = j(N*T,ncolumns,0);
DO i=1 TO N;
  MatrixCreateNew3[(i-1)*T+1:i*T, 1:ncolumns] = MatrixCreateNew2[1:T, 1:ncolumns];
END;



/*****************************************************************************************/
/*****************************************************************************************/
/*****************************************************************************************/
/*****************************************************************************************/


/*Keep backwards conditions*/
typeMtxfut=typeMtx3;
typeMtxFut[,loweridx]=0;

/*CHECKING HOW TO CHANGE ORDER OF COVARIATES FOR OUTCOMES AS PREDICTORS*/

CREATE TypeMtxFut from TypeMtxFut;
APPEND from TypeMtxFut;
CLOSE TypeMtxFut;



*Create the combination TypeMtx, of maximum size to be pared down later, excluding the intercept;
TypeMtxCombo = j((Np*T),T*T,0); 

TypeMtxCombo2 = j((Np*T),T*T,0);
*Loop across the different types, starting with current;
DO i = 0 TO (T-1);
	*Identify the indices for the appropriate setting;
	Idx = loc(rt-ct = i);
	IdxDif = setdif(1:T*T, Idx);/*THIS RETURNS ALL THOSE INDICES OF THE MATRIX FOR WHICH RT-CT IS DIFFERENT FROM I*/
	temp = typeMtxNB[2:(Np+1),];
	temp[,IdxDif] = 0;
	TypeMtxCombo[(i*Np+1):((i+1)*Np),] = temp;

	*Shift the values to accomodate shifted Y-X relationship;
	if(i>0) then DO; 
		dummyZero = j(Np, i, 0) || temp;
		TypeMtxCombo2[(i*Np+1):((i+1)*Np),] = dummyZero[,1:(T*T)];
	END;
	
	ELSE TypeMtxCombo2[(i*Np+1):((i+1)*Np),] = temp;

END; 



/***************************************************************************************************************/
/********** WILL HAVE TO CREATE SOMETHING SIMILAR TO TYPEMTXCOMBO 2 FOR EACH TIME-DEPENDENT COVARIATE***********/
/**************THIS USING TYPEMTXFUT****************************************************************************/
/***************************************************************************************************************/


*Identify predictors or relationships with no valid moment conditions;

neq = TypeMtxCombo2[,+]; /*ADDING UP NUMBER OF 1'S FOR EACH ROW OF TYPEMTXCOMBO2*/




*Print a note about omitted moment conditions;
if min(neq)=0 then do;
	zeroPred = loc(neq=0);
	Note = "There are no valid moment conditions for " +char(ncol(zeroPred))+" covariate relationship(s).";
	print Note[label="Moment Condition Notes"];
	print "These covariate relationships will be omitted in the analysis.";
end;
else print "All covariate relationships will be evaluated."[label="Moment Condition Notes"];


/*FOR THOSE ROWS RELATED TO LAGS OF THE COVARIATE WITH NO VALID MOMENTS WE DELETE THAT ROW*/
/*THIS BECAUSE WE WILL NOT BE USING THOSE LAGS IN OUR MODELS*/
keepPred = loc(neq>0);

TypeMtxCombo3 = TypeMtxCombo2[keepPred,]; /*WILL ONLY KEEP ROWS OF COVARIATES WITH VALID MOMENTS*/


*Append the intercept conditions to the matrix;
*Append the time independent conditions to the matrix, if available;
%if &predVarTI. ^=. %then %do;
/*	TypeMtxCombo3 = TypeMtxNB[1,] // TypeMtxIND // TypeMtxCombo3; *Tested intercept;*/
	TypeMtxCombo3 = shape(I(T),1) // TypeMtxIND // TypeMtxCombo3; *Type I intercept;
	
%end;
%else %do;
/*	TypeMtxCombo3 = TypeMtxNB[1,] // TypeMtxCombo3; *Tested intercept;*/
	TypeMtxCombo3 = shape(I(T),1) // TypeMtxCombo3; *Type I intercept;
%end;

*If there are time independent variables, they are represented after the intercept, before the time dependent;
CREATE TypeMtxCombo3 from TypeMtxCombo3;
APPEND from TypeMtxCombo3;
CLOSE TypeMtxCombo3;
print "Each row of TypeMtx is the shifted type vector for each of the predictors, by individual test";
print TypeMtxCombo3;

*Create a printable form of the type matrix;
Type2[,upperIdx] = 0;
Type4out = j(T,T*(Np),.);
DO i=1 to Np;
Type4out[,T*(i-1)+1:T*i] = shape(Type2[i,],T,T);
END;

CREATE Type4out from Type4out;
APPEND from Type4out;
CLOSE Type4out;


*Create the modified data set;
USE mydatasorted;
read all var {&predVarTD.} into xnew;
read all var {&outVar.} into Ynew;
read all var {&idVar. &timeVar.} into othersnew;
%if &predVarTI. ^=. %then %do;
	read all VAR {&idVar.  &timeVar. &outVar. &predVarTI.} INTO otherVars;
%end;
%else %do;
	read all VAR {&idVar.  &timeVar. &outVar.} INTO otherVars;
%end;
CLOSE mydatasorted;


*Creating the new  variables;
X2 = j(N*T,Np*(T-1),0);/*X2 has dimensions [total # obs, #of lagged predictors multiplied by number of lags, all values start at 0*/
DO i=1 TO T-1;
	X2[(i*N+1):N*T, ((i-1)*Np+1):i*Np] = Xnew[1:(T-i)*N,];
	/*from row N+1 to row total number of subjects, from column 1 to Np*/ 
    /*creates matrix with first Np columns containing lag 1 measurements, second Np columns containsing lag2 measurements*/

/*X2=[00|00  first two columns represent lag-1 measurments for both time-dependent covariates measured 3 times
      00|00  
	--------
	  11|00  Second two columns represent lag-2 measurements for both time-dependent covariates measured 3 times
	  11|00
	--------
	  22|11
	  22|11]  */
END;


X2 = Xnew || X2;


*Remove the predictors with no valid moments;
X3 = X2[,keepPred];


*Add in the ID, outcome and time (and possible time independent variables;
X3 = otherVars || X3;


*Create the adjusted variable names;
predNames = t(repeat(t({&predVarTD.}),T));
lagName = j(1,Np*T,0);
DO i=1 to (T-1);
	lagName[1,(i*Np+1):(i+1)*Np]=j(1,Np,i);
END;
varNames = catx("_", predNames, char(lagName));
varnamesTD=varnames;

*Retain only the variable names for those with valid moment conditions;
varNames = varNames[,keepPred];
%if &predVarTI. ^=. %then %do;
	varNames = {&predVarTI.} ||  varNames;
%end;

*Add in the names for the ID, outcome and time;
varNames2 = {&idVar.} || {&timeVar.} || {&outVar.} ||  varNames;

*Sort the data;
call sort(X3 , {1 2});

NcolX3=ncol(X3);

X3[1:N*T,4+NpTI2:NcolX3]=X3[1:N*T,4+NpTI2:NcolX3]#MatrixCreateNew3;


CREATE Mydata3 from X3[c=varNames2];
APPEND from X3;
close Mydata3;

%end;

QUIT;
%mend partitionedDataMatrix; *End macro code;




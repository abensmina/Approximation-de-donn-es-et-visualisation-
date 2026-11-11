	 PROGRAM FVSURFC1
      write(6,*)' '
      write(6,*)' '
      write(6,*)' '
      write(6,*)' '
      write(6,*)' '
      write(6,*)' '
      write(6,*)'*****************************************'
      write(6,*)'******  Programme SPLINE C1 usuel *******'
      write(6,*)'******  sans prise en compte      *******'
      write(6,*)'******  des fortes variations     *******'
      write(6,*)'*****************************************'
            write(6,*)' '
      write(6,*)' '
      write(6,*)' '
      write(6,*)' '

      CALL FVSURF1

      STOP
      END

c     *********************************************
c     *********************************************

      SUBROUTINE FVSURF1

      DOUBLE PRECISION X0,Y0,X1,Y1,X2,Y2,X3,Y3
      DOUBLE PRECISION AX,AY
      DOUBLE PRECISION PAHX,PAHY
      DOUBLE PRECISION FX(16,4),FY(16,4),DFX(16,4),DFY(16,4)
      DOUBLE PRECISION R,B,REF1(16,16),R1,EPS
      DOUBLE PRECISION DAFON,COTT
      DOUBLE PRECISION ALPHA,BETA

      CHARACTER*80 NOMFIC3
      CHARACTER*80 NOMFIC6
      CHARACTER*1 REP
      CHARACTER*1 REP1

      PARAMETER(MXNP=41000)
      PARAMETER(MXDIVX=30,MXDIVY=30)
      PARAMETER(MXIPAX=MXDIVX+1,MXIPAY=MXDIVY+1)
      PARAMETER(MXIDIM=4*MXIPAX*MXIPAY)

      DIMENSION AX(MXIPAX),AY(MXIPAY)
      DIMENSION R(MXIDIM,MXIDIM)
      DIMENSION B(MXIDIM)
      DIMENSION R1(MXIDIM,MXIDIM)
      DIMENSION DAFON(3,MXNP),COTT(MXNP)

      DIMENSION INC(MXDIVX,MXDIVY)
      DIMENSION INF(MXIPAX,MXIPAY,4)
      DIMENSION INDI(16,MXDIVX*MXDIVY)

      COMMON/A1/X0,Y0,X1,Y1,X2,Y2,X3,Y3
      COMMON/B1/IDIVX,IDIVY
      COMMON/C1/PAHX,PAHY
      COMMON/D1/IPAX,IPAY
      COMMON/E1/AX,AY
      COMMON/F1/INC,INDI
      COMMON/G1/INF
      COMMON/H1/IDIM
      COMMON/I1/REP1
      COMMON/K1/EPS

      DATA ALPHA/1.D+00/
       DATA BETA/0.000001/
c     ALPHA ecrase l'intervalle d' etude pour TRANSFO: 1 =(0,1) 
c     BETA ecrete les valeurs calculees



      OPEN (UNIT=1,FILE='s1don')
      OPEN (UNIT=2,FILE='c1ref')






1     CONTINUE
      WRITE(6,*) 'AJUSTEMENT DE DONNEES PONCTUELLES ?(o/n) :'

      READ(5,7) REP1

      IF(REP1.EQ.'n') GO TO 2


      IF(REP1.EQ.'o') THEN

        CALL DOMAINE

        CALL GRILLE

        READ(1,*) ((FX(L,I),I=1,4),L=1,16)
        READ(1,*) ((FY(L,I),I=1,4),L=1,16)
        READ(1,*) ((DFX(L,I),I=1,4),L=1,16)
        READ(1,*) ((DFY(L,I),I=1,4),L=1,16)
        CALL CFDBN(FX,FY,DFX,DFY)

        CALL NUMER

        CALL INITIAL(R,B,R1)

        IDIM6=IDIM

        WRITE(6,*) ' '
        WRITE(6,*) '* renseignements sur les donnes ponctuelles  ajust
     &er *'
        WRITE(6,*) ' '


          WRITE(6,*) '     NOM DU FICHIER CONTENANT LES POINTS DE DONNEE
     &S :'
          READ(5,7) NOMFIC3
          OPEN(UNIT=8,FILE=NOMFIC3)

        CALL AJDONP(FX,FY,R,B,NPAFON,DAFON,COTT,NJ)



      ELSE
        GO TO 1
      ENDIF




2     CONTINUE
      WRITE(6,*) ' '



7     FORMAT(A)




c                 ***********************



      READ(2,*) ((REF1(I,J),J=1,I),I=1,16)

      CALL REG(REF1,R1)

      CALL EPSI

      DO 10 K=1,IDIM
      DO 10 L=1,IDIM
        R(K,L)=R(K,L)+EPS*R1(K,L)
10    CONTINUE

       DO 11 I=1,IDIM6
         DO 12 J=1,I
            R(I,J)=R(J,I)
12       CONTINUE
11    CONTINUE


      CALL GAUSS(R,B,IDIM6)




c                 ***********************



      WRITE(6,*) ' '
      WRITE(6,*) '* estimation de la surface calcule *'
      WRITE(6,*) ' '

3     CONTINUE
      WRITE(6,*) 'ESTIMATION SUR UNE GRILLE REGULIERE ?(o/n) :'

      READ(5,7) REP

      IF(REP.EQ.'n') THEN
        GO TO 4
      ENDIF



      IF(REP.EQ.'o') THEN



        WRITE(6,*) '     NOM DU FICHIER RESULTAT :'
        READ(5,7) NOMFIC6
        OPEN(UNIT=11,FILE=NOMFIC6)


        CALL EVGRI(B,FX,FY,DFX,DFY,NPAFON,DAFON,COTT,NJ)

      ELSE
        GO TO 3
      ENDIF




4     CONTINUE

      CLOSE(UNIT=8)
      CLOSE(UNIT=11)



c                 ***********************


      CLOSE(UNIT=1)
      CLOSE(UNIT=2)

      WRITE(6,*) ' '
      WRITE(6,*) ' '
      WRITE(6,*) '                     *** fin d''excution du programme
     & ***'


      RETURN
      END

c     *********************************************
c     *********************************************

      SUBROUTINE DOMAINE

      DOUBLE PRECISION X0,Y0,X1,Y1,X2,Y2,X3,Y3
      DOUBLE PRECISION XMIN,XMAX,YMIN,YMAX

      COMMON/A1/X0,Y0,X1,Y1,X2,Y2,X3,Y3


      WRITE(6,*) ' '
      WRITE(6,*) '* renseignements concernant les limites de l''tude *'

      WRITE(6,*) ' '
      WRITE(6,*) 'COORDONNEES DU DOMAINE D''ETUDE :'
      WRITE(6,*) '     XMIN = '
      READ(5,*) XMIN
      WRITE(6,*) '     YMIN = '
      READ(5,*) YMIN
      WRITE(6,*) '     XMAX = '
      READ(5,*) XMAX
      WRITE(6,*) '     YMAX = '
      READ(5,*) YMAX


      X0=XMIN
      X1=XMAX
      X2=XMAX
      X3=XMIN
      Y0=YMIN
      Y1=YMIN
      Y2=YMAX
      Y3=YMAX



      RETURN
      END

c     *********************************************
c     *********************************************

      SUBROUTINE GRILLE

      DOUBLE PRECISION X0,Y0,X1,Y1,X2,Y2,X3,Y3
      DOUBLE PRECISION PAHX,PAHY
      DOUBLE PRECISION AX,AY

      PARAMETER(MXDIVX=30,MXDIVY=30)
      PARAMETER(MXIPAX=MXDIVX+1,MXIPAY=MXDIVY+1)

      DIMENSION AX(MXIPAX),AY(MXIPAY)

      COMMON/A1/X0,Y0,X1,Y1,X2,Y2,X3,Y3
      COMMON/B1/IDIVX,IDIVY
      COMMON/C1/PAHX,PAHY
      COMMON/D1/IPAX,IPAY
      COMMON/E1/AX,AY




      WRITE(6,*) ' '
      WRITE(6,*) '* paramtres de la grille de calcul *'

      WRITE(6,*) ' '
      WRITE(6,*) 'NOMBRE D''INTERVALLES DE DISCRETISATION EN X :'
      READ(5,*) IDIVX
      WRITE(6,*) 'NOMBRE D''INTERVALLES DE DISCRETISATION EN Y :'
      READ(5,*) IDIVY



      PAHX=(X1-X0)/DFLOAT(IDIVX)
      PAHY=(Y3-Y0)/DFLOAT(IDIVY)
      IPAX=IDIVX+1
      IPAY=IDIVY+1



      DO 1 I=1,IPAX
         AX(I)=X0+(DFLOAT(I-1))*PAHX
1     CONTINUE
      DO 2 J=1,IPAY
         AY(J)=Y0+(DFLOAT(J-1))*PAHY
2     CONTINUE



      RETURN
      END

c     *********************************************
c     *********************************************

      SUBROUTINE CFDBN(FX,FY,DFX,DFY)

      DOUBLE PRECISION FX(16,4),FY(16,4),DFX(16,4),DFY(16,4)
      DOUBLE PRECISION PAHX,PAHY

      COMMON/C1/PAHX,PAHY




      DO 1 LL=1,13,4
        DO 2 II=1,4
          FX(LL,II)=FX(LL,II)
          FY(LL,II)=FY(LL,II)
          DFX(LL,II)=(DFX(LL,II))/PAHX
          DFY(LL,II)=(DFY(LL,II))/PAHY
2       CONTINUE
1     CONTINUE

      DO 3 LL=2,14,4
        DO 4 II=1,4
          FX(LL,II)=PAHX*FX(LL,II)
          FY(LL,II)=FY(LL,II)
          DFY(LL,II)=(DFY(LL,II))/PAHY
4       CONTINUE
3     CONTINUE

      DO 5 LL=3,15,4
        DO 6 II=1,4
          FX(LL,II)=FX(LL,II)
          FY(LL,II)=PAHY*FY(LL,II)
          DFX(LL,II)=(DFX(LL,II))/PAHX
6       CONTINUE
5     CONTINUE


      DO 7 LL=4,16,4
        DO 8 II=1,4
          FX(LL,II)=PAHX*FX(LL,II)
          FY(LL,II)=PAHY*FY(LL,II)
8       CONTINUE
7     CONTINUE

      RETURN
      END

c     ********************************************
c     ********************************************

      SUBROUTINE  NUMER

      PARAMETER(MXDIVX=30,MXDIVY=30)
      PARAMETER(MXIPAX=MXDIVX+1,MXIPAY=MXDIVY+1)

      DIMENSION INC(MXDIVX,MXDIVY)
      DIMENSION INF(MXIPAX,MXIPAY,4)
      DIMENSION INDI(16,MXDIVX*MXDIVY)

      COMMON/B1/IDIVX,IDIVY
      COMMON/D1/IPAX,IPAY
      COMMON/F1/INC,INDI
      COMMON/G1/INF
      COMMON/H1/IDIM




      DO 1 M=1,IDIVY
        DO 2 L=1,IDIVX
           INC(L,M)=L+(M-1)*IDIVX
2       CONTINUE
1     CONTINUE

      DO 3 J=1,IPAY
        DO 4 I=1,IPAX
          DO 5 L=1,4
            INF(I,J,L)=(J-1)*4*IPAX+4*(I-1)+L
5         CONTINUE
4       CONTINUE
3     CONTINUE

      IDIM=INF(IPAX,IPAY,4)


      DO 6 J=1,IDIVY
         DO 7 I=1,IDIVX
            KI=INC(I,J)
            DO 8 L=1,4
              INDI(L,KI)=INF(I,J,L)
8           CONTINUE
            DO 9 L=5,8
              INDI(L,KI)=INF(I+1,J,L-4)
9           CONTINUE
            DO 10 L=9,12
              INDI(L,KI)=INF(I+1,J+1,L-8)
10          CONTINUE
            DO 11 L=13,16
              INDI(L,KI)=INF(I,J+1,L-12)
11          CONTINUE
7        CONTINUE
6     CONTINUE



      RETURN
      END

c     *********************************************
c     *********************************************

      SUBROUTINE INITIAL(R,B,R1)

      DOUBLE PRECISION R,B,R1

      PARAMETER(MXNP=41000)
      PARAMETER(MXDIVX=30,MXDIVY=30)
      PARAMETER(MXIPAX=MXDIVX+1,MXIPAY=MXDIVY+1)
      PARAMETER(MXIDIM=4*MXIPAX*MXIPAY)


      DIMENSION R(MXIDIM,MXIDIM)
      DIMENSION B(MXIDIM)
      DIMENSION R1(MXIDIM,MXIDIM)

      COMMON/H1/IDIM




      DO 1 K=1,MXIDIM
        DO 2 L=1,MXIDIM
          R(K,L)=0.D+00
2       CONTINUE
1     CONTINUE


      DO 3 K=1,MXIDIM
        B(K)=0.D+00
3     CONTINUE


      DO 4 K=1,IDIM
        DO 5 L=1,IDIM
          R1(K,L)=0.D+00
5       CONTINUE
4     CONTINUE



      RETURN
      END

c     *********************************************
c     *********************************************

      SUBROUTINE REG(REF1,R1)

      DOUBLE PRECISION PAHX,PAHY
      DOUBLE PRECISION REF1(16,16),R1
      DOUBLE PRECISION C(16)

      PARAMETER(MXDIVX=30,MXDIVY=30)
      PARAMETER(MXIPAX=MXDIVX+1,MXIPAY=MXDIVY+1)
      PARAMETER(MXIDIM=4*MXIPAX*MXIPAY)

      DIMENSION R1(MXIDIM,MXIDIM)
      DIMENSION INC(MXDIVX,MXDIVY)
      DIMENSION INDI(16,MXDIVX*MXDIVY)


      COMMON/B1/IDIVX,IDIVY
      COMMON/C1/PAHX,PAHY
      COMMON/F1/INC,INDI





      DO 1 L=1,13,4
1       C(L)=1.D+00
      DO 2 L=2,14,4
2       C(L)=PAHX
      DO 3 L=3,15,4
3       C(L)=PAHY
      DO 4 L=4,16,4
4       C(L)=PAHX*PAHY

      DO 5 I=1,16
        DO 6 J=1,16
          REF1(I,J)=REF1(J,I)
6       CONTINUE
5     CONTINUE


      DO 7 J=1,IDIVY
         DO 8 I=1,IDIVX
            KI=INC(I,J)
            DO 9 L=1,16
               DO 10 M=1,16
                  R1(INDI(L,KI),INDI(M,KI))=R1(INDI(L,KI),INDI(M,KI))
     &                              +(C(L)*C(M)*REF1(L,M))*PAHX*PAHY
10             CONTINUE
9           CONTINUE
8        CONTINUE
7     CONTINUE



      RETURN
      END

c     *********************************************
c     *********************************************

      SUBROUTINE EPSI

      DOUBLE PRECISION EPS
      CHARACTER*1 REP
      CHARACTER*1 REP1

      COMMON/I1/REP1
      COMMON/K1/EPS



      IF(REP1.EQ.'o') THEN
        EPS=1.D-06
      ENDIF

1     CONTINUE
        WRITE(6,*) ' '
        WRITE(6,*) '* choix particulier du paramtre de lissage *'
        WRITE(6,*) ' '
      WRITE(6,*) 'VOULEZ-VOUS UNE VALEUR PARTICULIERE DU PARAMETRE DE LI
     &SSAGE ?(o/n) :'
      READ(5,2) REP
      IF(REP.EQ.'n') THEN
        GO TO 3
      ENDIF
      IF(REP.EQ.'o') THEN
        WRITE(6,*) '     VALEUR DU PARAMETRE DE LISSAGE :'
        READ(5,*) EPS
      ELSE
        GO TO 1
      ENDIF

2     FORMAT(A)
3     CONTINUE




      RETURN
      END


c     *********************************************
c     *********************************************

      SUBROUTINE AJDONP(FX,FY,R,B,NPAFON,DAFON,COTT,NJ)

      DOUBLE PRECISION DAFON
      DOUBLE PRECISION COTT
      DOUBLE PRECISION FX(16,4),FY(16,4)
      DOUBLE PRECISION R,B
      DOUBLE PRECISION VV


      PARAMETER(MXNP=41000)
      PARAMETER(MXDIVX=30,MXDIVY=30)
      PARAMETER(MXIPAX=MXDIVX+1,MXIPAY=MXDIVY+1)
      PARAMETER(MXIDIM=4*MXIPAX*MXIPAY)

      DIMENSION DAFON(3,MXNP)
      DIMENSION COTT(MXNP)
      DIMENSION R(MXIDIM,MXIDIM)
      DIMENSION B(MXIDIM)
      DIMENSION VV(MXNP,MXIDIM)

      COMMON/H1/IDIM





      CALL PARAMP(NPAFON,DAFON,COTT,NJ)

      CALL VALEUR(VV,FX,FY,NPAFON,DAFON)
      
      write(6,*)'npafon',npafon


      DO 1 J=1,NPAFON
        DO 2 IK=1,IDIM
          B(IK)=B(IK)+COTT(J)*VV(J,IK)
          DO 3 IL=1,IDIM
            R(IK,IL)=R(IK,IL)+VV(J,IK)*VV(J,IL)
3         CONTINUE
2       CONTINUE
1     CONTINUE




      RETURN
      END


c     *********************************************
c     *********************************************

      SUBROUTINE VALEUR(VVF,XF,YF,NP,XX)

      DOUBLE PRECISION PAHX,PAHY
      DOUBLE PRECISION XX,VVF
      DOUBLE PRECISION D1X,D2X,D1Y,D2Y,XXX1,XXX2
      DOUBLE PRECISION AX,AY
      DOUBLE PRECISION P(4),Q(4)
      DOUBLE PRECISION XF(16,4),YF(16,4)
      DOUBLE PRECISION S1,S2

      PARAMETER(MXNP=41000)
      PARAMETER(MXDIVX=30,MXDIVY=30)
      PARAMETER(MXIPAX=MXDIVX+1,MXIPAY=MXDIVY+1)
      PARAMETER(MXIDIM=4*MXIPAX*MXIPAY)

      DIMENSION XX(3,MXNP),VVF(MXNP,MXIDIM)
      DIMENSION AX(MXIPAX),AY(MXIPAY)
      DIMENSION INC(MXDIVX,MXDIVY)
      DIMENSION INDI(16,MXDIVX*MXDIVY)

      COMMON/B1/IDIVX,IDIVY
      COMMON/C1/PAHX,PAHY
      COMMON/E1/AX,AY
      COMMON/F1/INC,INDI
      COMMON/H1/IDIM





      DO 1 J=1,NP
        DO 2 K=1,IDIM
          VVF(J,K)=0.D+00
2       CONTINUE
1     CONTINUE



      DO 3 J=1,NP
        DO 4 M=1,IDIVY
          DO 5 L=1,IDIVX
            KI=INC(L,M)
            D1X=XX(1,J)-AX(L)
            D2X=XX(1,J)-AX(L+1)
            D1Y=XX(2,J)-AY(M)
            D2Y=XX(2,J)-AY(M+1)
            XXX1=D1X*D2X
            XXX2=D1Y*D2Y

            IF ((XXX1.LE.(0.)).AND.(XXX2.LE.(0.))) THEN
              DO 6 LL=1,16
                DO 7 II=1,4
                  P(II)=XF(LL,II)
                  Q(II)=YF(LL,II)
7               CONTINUE
                CALL VAL(P,XX(1,J),AX(L),PAHX,S1)

                CALL VAL(Q,XX(2,J),AY(M),PAHY,S2)

                VVF(J,INDI(LL,KI))=VVF(J,INDI(LL,KI))+S1*S2
6             CONTINUE
            GO TO 8
            ENDIF

5         CONTINUE
4       CONTINUE
8       CONTINUE
3     CONTINUE

      RETURN
      END

c     *********************************************
c     *********************************************

      SUBROUTINE PARAMP(NPAFON,DAFON,COTT,NJ)

      DOUBLE PRECISION DAFON
      DOUBLE PRECISION DDAFON
      DOUBLE PRECISION COTT

      PARAMETER(MXNP=41000)

      DIMENSION DAFON(3,MXNP)
      DIMENSION DDAFON(3,MXNP)
      DIMENSION COTT(MXNP)


      READ(8,*) NPAFON
      DO 1 J=1,NPAFON
        READ(8,*) DDAFON(1,J),DDAFON(2,J),DDAFON(3,J)
         DAFON(1,J)=DDAFON(1,J)
         DAFON(2,J)=DDAFON(2,J)
         DAFON(3,J)=DDAFON(3,J)
        COTT(J)=DAFON(3,J)
1     CONTINUE




      RETURN
      END

c     *********************************************
c     *********************************************

      SUBROUTINE GAUSS(R,B,IDIM6)

      DOUBLE PRECISION R,B
      DOUBLE PRECISION MAX,CC,T,S

      PARAMETER(MXNP=41000)
      PARAMETER(MXDIVX=30,MXDIVY=30)
      PARAMETER(MXIPAX=MXDIVX+1,MXIPAY=MXDIVY+1)
      PARAMETER(MXIDIM=4*MXIPAX*MXIPAY)

      DIMENSION R(MXIDIM,MXIDIM)
      DIMENSION B(MXIDIM)





      DO 1 K=1,IDIM6-1
        MAX=DABS(R(K,K))
        INDICE=K
        DO 2 L=K+1,IDIM6
          IF(DABS(R(L,K)).GT.MAX) THEN
            MAX=DABS(R(L,K))
            INDICE=L
          ENDIF
2       CONTINUE
        IF(INDICE.NE.K) THEN
          DO 3 L=K,IDIM6
            CC=R(K,L)
            R(K,L)=R(INDICE,L)
            R(INDICE,L)=CC
3         CONTINUE
          CC=B(K)
          B(K)=B(INDICE)
          B(INDICE)=CC
        ENDIF
        DO 4 L=K+1,IDIM6
          T=(R(L,K))/(R(K,K))
          DO 5 J=K,IDIM6
            R(L,J)=R(L,J)-(T*R(K,J))
5         CONTINUE
          B(L)=B(L)-(T*B(K))
4       CONTINUE
1     CONTINUE

      B(IDIM6)=(B(IDIM6))/(R(IDIM6,IDIM6))
      DO 6 I=IDIM6-1,1,-1
        S=B(I)
        DO 7 J=I+1,IDIM6
          S=S-B(J)*R(I,J)
7       CONTINUE
        B(I)=S/R(I,I)
6     CONTINUE




      RETURN
      END

c     *********************************************
c     *********************************************

      SUBROUTINE EVA(U,V,V1,XF,YF,BB)

      DOUBLE PRECISION V1,U,V
      DOUBLE PRECISION PAHX,PAHY
      DOUBLE PRECISION AX,AY
      DOUBLE PRECISION D1X,D2X,D1Y,D2Y,XXX1,XXX2
      DOUBLE PRECISION XF(16,4),YF(16,4)
      DOUBLE PRECISION P(4),Q(4),S1,S2
      DOUBLE PRECISION BB

      PARAMETER(MXNP=41000)
      PARAMETER(MXDIVX=30,MXDIVY=30)
      PARAMETER(MXIPAX=MXDIVX+1,MXIPAY=MXDIVY+1)
      PARAMETER(MXIDIM=4*MXIPAX*MXIPAY)

      DIMENSION AX(MXIPAX),AY(MXIPAY)
      DIMENSION BB(MXIDIM)
      DIMENSION INC(MXDIVX,MXDIVY)
      DIMENSION INDI(16,MXDIVX*MXDIVY)

      COMMON/B1/IDIVX,IDIVY
      COMMON/C1/PAHX,PAHY
      COMMON/E1/AX,AY
      COMMON/F1/INC,INDI





      V1=0.D+00

      DO 1 M=1,IDIVY
        DO 2 L=1,IDIVX

          KI=INC(L,M)
          D1X=U-AX(L)
          D2X=U-AX(L+1)
          D1Y=V-AY(M)
          D2Y=V-AY(M+1)
          XXX1=D1X*D2X
          XXX2=D1Y*D2Y

          IF ((XXX1.LE.(0.)).AND.(XXX2.LE.(0.))) THEN
            DO 3 LL=1,16
              DO 4 II=1,4
                P(II)=XF(LL,II)
                Q(II)=YF(LL,II)
4             CONTINUE

              CALL VAL(P,U,AX(L),PAHX,S1)
              CALL VAL(Q,V,AY(M),PAHY,S2)

              V1=V1+BB(INDI(LL,KI))*S1*S2

3           CONTINUE
            GO TO 5
          ENDIF

2       CONTINUE
1     CONTINUE
5     CONTINUE





      RETURN
      END

c     *********************************************
c     *********************************************

      SUBROUTINE VAL(PP,D,PA,PAH,S)

      DOUBLE PRECISION PP(4),D,PA,PAH,S,U


      U=(D-PA)/PAH
      S=PP(1)+U*(PP(2)+U*(PP(3)+U*PP(4)))



      RETURN
      END

c     *********************************************
c     *********************************************

      SUBROUTINE EVGRI(B,FX,FY,DFX,DFY,NPAFON,DAFON,COTT,NJ)


      DOUBLE PRECISION X0,Y0,X1,Y1,X2,Y2,X3,Y3
      DOUBLE PRECISION PAHX,PAHY
      DOUBLE PRECISION AX,AY
      DOUBLE PRECISION FX(16,4),FY(16,4)
      DOUBLE PRECISION DFX(16,4),DFY(16,4)
      DOUBLE PRECISION B
      DOUBLE PRECISION PAEVX,PAEVY,EVX,EVY
      DOUBLE PRECISION U,V,VLIS
      DOUBLE PRECISION DAFON,COTT


      PARAMETER(MXNP=41000)
      PARAMETER(MXDIVX=30,MXDIVY=30)
      PARAMETER(MXNPEVX=1000,MXNPEVY=1000)
      PARAMETER(MXIPAX=MXDIVX+1,MXIPAY=MXDIVY+1)
      PARAMETER(MXIDIM=4*MXIPAX*MXIPAY)
	

      DIMENSION AX(MXIPAX),AY(MXIPAY)
      DIMENSION B(MXIDIM)
      DIMENSION EVX(MXNPEVX),EVY(MXNPEVY)
      DIMENSION VLIS(MXNPEVX)
      DIMENSION DAFON(3,MXNP)
      DIMENSION COTT(MXNP)

      DIMENSION INC(MXDIVX,MXDIVY)
      DIMENSION INDI(16,MXDIVX*MXDIVY)

      COMMON/A1/X0,Y0,X1,Y1,X2,Y2,X3,Y3
      COMMON/B1/IDIVX,IDIVY
      COMMON/C1/PAHX,PAHY
      COMMON/E1/AX,AY
      COMMON/F1/INC,INDI




      WRITE(6,*) '     NOMBRE DE POINTS D''EVALUATION EN X :'
      READ(5,*) NPEVX
      WRITE(6,*) '     NOMBRE DE POINTS D''EVALUATION EN Y :'
      READ(5,*) NPEVY


      WRITE(11,*) 'VARIABLES="X","Y","Z"'
      WRITE(11,*) 'ZONE I=',NPEVX,',J=',NPEVY,',F=POINT'

      PAEVX=(X1-X0)/DFLOAT(NPEVX-1)
      PAEVY=(Y3-Y0)/DFLOAT(NPEVY-1)



      DO 1 J=1,NPEVY
        DO 2 I=1,NPEVX
          EVX(I)=X0+(DFLOAT(I-1))*PAEVX
          EVY(J)=Y0+(DFLOAT(J-1))*PAEVY
          U=EVX(I)
          V=EVY(J)

          CALL EVA(U,V,VLIS(I),FX,FY,B)

          WRITE(11,*) U,V,(VLIS(I))




2       CONTINUE
1     CONTINUE


c3     FORMAT(3(F18.6))



      RETURN
      END

c     *********************************************
c     *********************************************
c     *********************************************
c     *********************************************

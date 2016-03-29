#include <Rcpp.h>
using namespace Rcpp;

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

//get the list element named str, or return NULL.
// SEXP getListElement(SEXP list, const char *str)
// {
//   SEXP elmt=R_NilValue,names=Rf_getAttrib(list,R_NamesSymbol);
//   int i;
//   for (i=0;i<length(list);i++)
//     if(strcmp(CHAR(STRING_ELT(names,i)),str)==0)
//     {
//       elmt=VECTOR_ELT(list,i);
//       break;
//     }
//     return elmt;
// }

//generate multivariate normal observations from a DAG.
// [[Rcpp::export]]
SEXP datGenC(SEXP node,SEXP maxdeg,SEXP ordex,SEXP ts,SEXP noi,SEXP nobs,SEXP ivn,SEXP fX,SEXP beta)
{
  int i,j,k,jx,s,xnode,xmaxdeg,*xordex,*xts,xnoi,xnobs,*xivn,*ind,nProtected=0,tmp1,tmp2,tmp3,pa;
  double *xfX,*xbeta,*rdata,pp;
  SEXP data;

  xnode=INTEGER(node)[0];
  xmaxdeg=INTEGER(maxdeg)[0];
  xordex=INTEGER(ordex);
  xts=INTEGER(ts);
  xnoi=INTEGER(noi)[0];
  xnobs=INTEGER(nobs)[0];
  xivn=INTEGER(ivn);
  xfX=REAL(fX);
  xbeta=REAL(beta);
  ind=(int *) Calloc(xnode,int);
  assert(ind);

  for (k=0;k<xnode;k++)
  {
    for (i=k*xmaxdeg;i<(k+1)*xmaxdeg;i++)
    {
      if(xordex[i]!=0)
      {
        ind[k]=1;
        break;
      }
    }
  }


  PROTECT(data=Rf_allocVector(REALSXP, xnode*xnobs*xnoi));
  rdata=REAL(data);
  ++nProtected;

  GetRNGstate();
  for (k=0;k<xnoi;k++)
  {
    tmp1=k*xnobs;
    tmp2=tmp1*xnode;
    for (i=0;i<xnobs;i++)
    {
      tmp3=i*xnode;
      for (j=0;j<xnode;j++)
      {
        jx=xts[j];
        if(jx==xivn[k])
        {
          rdata[tmp2+tmp3+jx-1]=xfX[tmp1+i];
        }
        else if(ind[jx-1]==0)
        {
          rdata[tmp2+tmp3+jx-1]=R::rnorm(0.0,1.0);
        }
        else
        {
          pp=0.0;
          for (s=0;s<xmaxdeg;s++)
          {
            pa=xordex[(jx-1)*xmaxdeg+s];
            if(pa!=0)
            {
              pp+=rdata[tmp2+tmp3+pa-1]*xbeta[(jx-1)*xmaxdeg+s];
            }
          }
          rdata[tmp2+tmp3+jx-1]=pp+R::rnorm(0.0,1.0);
        }
      }
    }
  }
  PutRNGstate();

  Free(ind);
  UNPROTECT(nProtected);
  return data;
}




//generate inprod and fprod.

// [[Rcpp::export]]
SEXP prodG(SEXP xnode,SEXP xlist,SEXP ylist,SEXP obslength)
{
  int i,j,k,l,node,node2,node3,nProtected=0,len,ind1,ind2,ind3,ind4;
  double *xtmp,*ytmp,*inprod_v,*fprod_v;
  node=INTEGER(xnode)[0];
  node2=node*node;
  node3=node*node*node;

  //creat the return list
  SEXP retList;
  PROTECT(retList=Rf_allocVector(VECSXP,2));
  ++nProtected;

  //creat and set the names attribute of the return list
  SEXP names;
  PROTECT(names=Rf_allocVector(STRSXP,2));
  ++nProtected;
  SET_STRING_ELT(names,0,Rf_mkChar("inprod"));
  SET_STRING_ELT(names,1,Rf_mkChar("fprod"));
  Rf_setAttrib(retList,R_NamesSymbol,names);

  //creat the elements of the return list
  SEXP inprod;
  PROTECT(inprod=Rf_allocVector(REALSXP,node2));
  inprod_v=REAL(inprod);
  ++nProtected;
  SEXP fprod;
  PROTECT(fprod=Rf_allocVector(REALSXP,node3));
  fprod_v=REAL(fprod);
  ++nProtected;

  ind1=0;
  for (i=0;i<node;i++)
  {
    ind2=0;
    xtmp=REAL(VECTOR_ELT(xlist,i));
    ytmp=REAL(VECTOR_ELT(ylist,i));
    len=INTEGER(obslength)[i];
    for (j=0;j<node;j++)
    {
      inprod_v[ind1]=0;
      for(k=0;k<len;k++)
      {
        inprod_v[ind1]+=xtmp[ind2]*ytmp[k];
        ind2++;
      }
      ind1++;
    }
  }

  for (i=0;i<node;i++)
  {
    xtmp=REAL(VECTOR_ELT(xlist,i));
    len=INTEGER(obslength)[i];
    ind1=i*node2;
    for (j=0;j<node;j++)
    {
      ind2=j*len;
      for (k=j;k<node;k++)
      {
        ind3=j*len;
        ind4=ind1+j*node+k;
        fprod_v[ind4]=0;
        for(l=0;l<len;l++)
        {
          fprod_v[ind4]+=xtmp[ind2]*xtmp[ind3];
          ind2++;
          ind3++;
        }
        fprod_v[ind1+k*node+j]=fprod_v[ind4];
      }
    }
  }

  //set the elements of the return list
  SET_VECTOR_ELT(retList,0,inprod);
  SET_VECTOR_ELT(retList,1,fprod);

  UNPROTECT(nProtected);
  return retList;
}

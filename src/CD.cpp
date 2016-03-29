#include <Rcpp.h>
using namespace Rcpp;

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


//return the maximum of two numbers.
double Max(double x,double y)
{
  if(x>y) return x;
  else return y;
}



//check whether adding a directed edge from node a to node b induces cycles in graph G.
int Cycle(int node,int *G, int a, int b)
{
  if(a==b)  return 1;
  int i,j,lo,nBot=0,nTop=0,SLeng=1,bCycle=0;
  int *color,*S;
  color=(int *) Calloc(node,int);
  S=(int *) Calloc(node,int);
  color[a-1]=1;
  S[0]=a;
  while(SLeng>0)
  {
    i=S[nBot];
    SLeng--;
    nBot++;
    for(j=0;j<node;j++)
    {
      lo=(i-1)*node;
      if(G[lo+j]==1)
      {
        if((j+1)==b)
        {
          bCycle=1;
          break;
        }
        else if(color[j]==0)
        {
          nTop++;
          S[nTop]=j+1;
          SLeng++;
          color[j]=1;
        }
      }
    }
    if(bCycle==1) break;
  }
  Free(color);
  Free(S);
  return bCycle;
}



//coordinate descent for a single lambda.
void CoordinateDescent(int node,double *beta,double lambda,double *inprod,double *fprod,int *eorder,double *rsq,double eps,int *iter,double *sigma,double *norm,int *obsleng,int *eorleng)
{
  int i,j,k,node2=node*node,rn,pn,rindex,pindex,active_length;
  int *G,*tG,*active,*eor;
  double b1,b2,cbeta1,cbeta2,dif1=0,dif2=0,rsq1=*rsq,rsq2=*rsq,mad;
  double denom1,denom2;
  double rnorm1,rnorm2,rsig1,rsig2,pnorm1,pnorm2,psig1,psig2,
  absum1,absum2,tmpSigSq1,tmpSigSq2,bda1,bda2,
  Delta1,Delta2,root1,root2,root3,tmpCoef1,tmpCoef2,tmpCoef3,dLoss1,dLoss2,dLoss3;
  G=(int *) Calloc(node2,int);
  tG=(int *) Calloc(node2,int);
  for (i=0;i<node2;i++)
  {
    if(beta[i]!=0)
      G[i]=1;
    else
      G[i]=0;
  }
  while(1)
  {
    *iter+=1;
    active=(int *) Calloc((node2-node)/2,int);
    active_length=0;
    mad=0.0;
    eor=eorder;
    for(i=0;i<*eorleng;i++)
    {
      rn=*eor++;
      pn=*eor++;
      pindex=(pn-1)*node+rn-1;
      rindex=(rn-1)*node+pn-1;
      for (j=0;j<node2;j++)
        tG[j]=G[j];
      tG[rindex]=tG[pindex]=0;
      b1=beta[rindex];
      b2=beta[pindex];

      if(Cycle(node,tG,rn,pn))
      {
        denom1=fprod[(rn-1)*node2+(pn-1)*(node+1)];
        rnorm2=norm[rn-1]+2*b1*inprod[rindex]+denom1*b1*b1;
        absum1=inprod[rindex]+denom1*b1;
        tmpSigSq1=rnorm2/obsleng[rn-1];
        bda1=0.5*obsleng[rn-1]*denom1/(lambda*lambda);
        Delta1=bda1*bda1-(denom1*rnorm2-absum1*absum1)/(lambda*lambda);
        if(absum1>lambda*tmpSigSq1)
        {
          rsig1=bda1-sqrt(Delta1);
          cbeta1=(absum1-lambda*rsig1)/denom1;
        }
        else if(absum1<(-lambda*tmpSigSq1))
        {
          rsig1=bda1-sqrt(Delta1);
          cbeta1=(absum1+lambda*rsig1)/denom1;
        }
        else if(Delta1>=0 && fabs(absum1)>lambda*bda1)
        {
          root1=bda1+sqrt(Delta1);
          root2=bda1-sqrt(Delta1);
          root3=tmpSigSq1;
          if(absum1>0)
          {
            tmpCoef1=(absum1-lambda*root1)/denom1;
            tmpCoef2=(absum1-lambda*root2)/denom1;
            tmpCoef3=0;
          } else
          {
            tmpCoef1=(absum1+lambda*root1)/denom1;
            tmpCoef2=(absum1+lambda*root2)/denom1;
            tmpCoef3=0;
          }
          dLoss1=0.5*obsleng[rn-1]*log(root1)+lambda*fabs(tmpCoef1);
          dLoss2=0.5*obsleng[rn-1]*log(root2)+lambda*fabs(tmpCoef2);
          dLoss3=0.5*obsleng[rn-1]*log(root3)+lambda*fabs(tmpCoef3);
          rsig1= (dLoss1<dLoss2) ? ((dLoss1<dLoss3) ? root1 : root3) : ((dLoss2<dLoss3) ? root2 : root3);
          cbeta1= (dLoss1<dLoss2) ? ((dLoss1<dLoss3) ? tmpCoef1 : tmpCoef3) : ((dLoss2<dLoss3) ? tmpCoef2 : tmpCoef3);
        }
        else
        {
          rsig1=tmpSigSq1;
          cbeta1=0;
        }
        dif1=cbeta1-b1;
        if(dif1!=0)
        {
          rnorm1=norm[rn-1]-2*dif1*inprod[rindex]+denom1*dif1*dif1;
          *rsq=*rsq+0.5*obsleng[rn-1]*log(rsig1/sigma[rn-1])+lambda*(fabs(cbeta1)-fabs(b1));
          beta[rindex]=cbeta1;
          norm[rn-1]=rnorm1;
          sigma[rn-1]=rsig1;
          for (k=(rn-1)*node;k<rn*node;k++)
          {
            inprod[k]-=fprod[(rn-1)*(node2-node)+(pn-1)*node+k]*dif1;
          }
        }
        if(cbeta1!=0)
        {
          G[rindex]=1;
          active[active_length]=i;
          active_length++;
        }
        else
          G[rindex]=0;
        mad=Max(fabs(dif1),mad);
      }

      else if(Cycle(node,tG,pn,rn))
      {
        denom2=fprod[(pn-1)*node2+(rn-1)*(node+1)];
        pnorm1=norm[pn-1]+2*b2*inprod[pindex]+denom2*b2*b2;
        absum2=inprod[pindex]+denom2*b2;
        tmpSigSq2=pnorm1/obsleng[pn-1];
        bda2=0.5*obsleng[pn-1]*denom2/(lambda*lambda);
        Delta2=bda2*bda2-(denom2*pnorm1-absum2*absum2)/(lambda*lambda);
        if(absum2>lambda*tmpSigSq2)
        {
          psig2=bda2-sqrt(Delta2);
          cbeta2=(absum2-lambda*psig2)/denom2;
        }
        else if(absum2<(-lambda*tmpSigSq2))
        {
          psig2=bda2-sqrt(Delta2);
          cbeta2=(absum2+lambda*psig2)/denom2;
        }
        else if(Delta2>=0 && fabs(absum2)>lambda*bda2)
        {
          root1=bda2+sqrt(Delta2);
          root2=bda2-sqrt(Delta2);
          root3=tmpSigSq2;
          if(absum2>0)
          {
            tmpCoef1=(absum2-lambda*root1)/denom2;
            tmpCoef2=(absum2-lambda*root2)/denom2;
            tmpCoef3=0;
          } else
          {
            tmpCoef1=(absum2+lambda*root1)/denom2;
            tmpCoef2=(absum2+lambda*root2)/denom2;
            tmpCoef3=0;
          }
          dLoss1=0.5*obsleng[pn-1]*log(root1)+lambda*fabs(tmpCoef1);
          dLoss2=0.5*obsleng[pn-1]*log(root2)+lambda*fabs(tmpCoef2);
          dLoss3=0.5*obsleng[pn-1]*log(root3)+lambda*fabs(tmpCoef3);
          psig2= (dLoss1<dLoss2) ? ((dLoss1<dLoss3) ? root1 : root3) : ((dLoss2<dLoss3) ? root2 : root3);
          cbeta2= (dLoss1<dLoss2) ? ((dLoss1<dLoss3) ? tmpCoef1 : tmpCoef3) : ((dLoss2<dLoss3) ? tmpCoef2 : tmpCoef3);
        }
        else
        {
          psig2=tmpSigSq2;
          cbeta2=0;
        }
        dif2=cbeta2-b2;
        if(dif2!=0)
        {
          pnorm2=norm[pn-1]-2*dif2*inprod[pindex]+denom2*dif2*dif2;
          *rsq=*rsq+0.5*obsleng[pn-1]*log(psig2/sigma[pn-1])+lambda*(fabs(cbeta2)-fabs(b2));
          beta[pindex]=cbeta2;
          norm[pn-1]=pnorm2;
          sigma[pn-1]=psig2;
          for (k=(pn-1)*node;k<pn*node;k++)
          {
            inprod[k]-=fprod[(pn-1)*(node2-node)+(rn-1)*node+k]*dif2;
          }
        }
        if(cbeta2!=0)
        {
          G[pindex]=1;
          active[active_length]=i;
          active_length++;
        }
        else
          G[pindex]=0;
        mad=Max(fabs(dif2),mad);
      }

      else
      {
        denom1=fprod[(rn-1)*node2+(pn-1)*(node+1)];
        denom2=fprod[(pn-1)*node2+(rn-1)*(node+1)];

        rnorm2=norm[rn-1]+2*b1*inprod[rindex]+denom1*b1*b1;
        absum1=inprod[rindex]+denom1*b1;
        tmpSigSq1=rnorm2/obsleng[rn-1];
        bda1=0.5*obsleng[rn-1]*denom1/(lambda*lambda);
        Delta1=bda1*bda1-(denom1*rnorm2-absum1*absum1)/(lambda*lambda);
        if(absum1>lambda*tmpSigSq1)
        {
          rsig1=bda1-sqrt(Delta1);
          cbeta1=(absum1-lambda*rsig1)/denom1;
        }
        else if(absum1<(-lambda*tmpSigSq1))
        {
          rsig1=bda1-sqrt(Delta1);
          cbeta1=(absum1+lambda*rsig1)/denom1;
        }
        else if(Delta1>=0 && fabs(absum1)>lambda*bda1)
        {
          root1=bda1+sqrt(Delta1);
          root2=bda1-sqrt(Delta1);
          root3=tmpSigSq1;
          if(absum1>0)
          {
            tmpCoef1=(absum1-lambda*root1)/denom1;
            tmpCoef2=(absum1-lambda*root2)/denom1;
            tmpCoef3=0;
          } else
          {
            tmpCoef1=(absum1+lambda*root1)/denom1;
            tmpCoef2=(absum1+lambda*root2)/denom1;
            tmpCoef3=0;
          }
          dLoss1=0.5*obsleng[rn-1]*log(root1)+lambda*fabs(tmpCoef1);
          dLoss2=0.5*obsleng[rn-1]*log(root2)+lambda*fabs(tmpCoef2);
          dLoss3=0.5*obsleng[rn-1]*log(root3)+lambda*fabs(tmpCoef3);
          rsig1= (dLoss1<dLoss2) ? ((dLoss1<dLoss3) ? root1 : root3) : ((dLoss2<dLoss3) ? root2 : root3);
          cbeta1= (dLoss1<dLoss2) ? ((dLoss1<dLoss3) ? tmpCoef1 : tmpCoef3) : ((dLoss2<dLoss3) ? tmpCoef2 : tmpCoef3);
        }
        else
        {
          rsig1=tmpSigSq1;
          cbeta1=0;
        }
        dif1=cbeta1-b1;
        rnorm1=norm[rn-1]-2*dif1*inprod[rindex]+denom1*dif1*dif1;

        pnorm1=norm[pn-1]+2*b2*inprod[pindex]+denom2*b2*b2;
        absum2=inprod[pindex]+denom2*b2;
        tmpSigSq2=pnorm1/obsleng[pn-1];
        bda2=0.5*obsleng[pn-1]*denom2/(lambda*lambda);
        Delta2=bda2*bda2-(denom2*pnorm1-absum2*absum2)/(lambda*lambda);
        if(absum2>lambda*tmpSigSq2)
        {
          psig2=bda2-sqrt(Delta2);
          cbeta2=(absum2-lambda*psig2)/denom2;
        }
        else if(absum2<(-lambda*tmpSigSq2))
        {
          psig2=bda2-sqrt(Delta2);
          cbeta2=(absum2+lambda*psig2)/denom2;
        }
        else if(Delta2>=0 && fabs(absum2)>lambda*bda2)
        {
          root1=bda2+sqrt(Delta2);
          root2=bda2-sqrt(Delta2);
          root3=tmpSigSq2;
          if(absum2>0)
          {
            tmpCoef1=(absum2-lambda*root1)/denom2;
            tmpCoef2=(absum2-lambda*root2)/denom2;
            tmpCoef3=0;
          } else
          {
            tmpCoef1=(absum2+lambda*root1)/denom2;
            tmpCoef2=(absum2+lambda*root2)/denom2;
            tmpCoef3=0;
          }
          dLoss1=0.5*obsleng[pn-1]*log(root1)+lambda*fabs(tmpCoef1);
          dLoss2=0.5*obsleng[pn-1]*log(root2)+lambda*fabs(tmpCoef2);
          dLoss3=0.5*obsleng[pn-1]*log(root3)+lambda*fabs(tmpCoef3);
          psig2= (dLoss1<dLoss2) ? ((dLoss1<dLoss3) ? root1 : root3) : ((dLoss2<dLoss3) ? root2 : root3);
          cbeta2= (dLoss1<dLoss2) ? ((dLoss1<dLoss3) ? tmpCoef1 : tmpCoef3) : ((dLoss2<dLoss3) ? tmpCoef2 : tmpCoef3);
        }
        else
        {
          psig2=tmpSigSq2;
          cbeta2=0;
        }
        dif2=cbeta2-b2;
        pnorm2=norm[pn-1]-2*dif2*inprod[pindex]+denom2*dif2*dif2;

        psig1=pnorm1/obsleng[pn-1];
        rsig2=rnorm2/obsleng[rn-1];

        rsq1=*rsq+0.5*obsleng[rn-1]*log(rsig1/sigma[rn-1])+lambda*(fabs(cbeta1)-fabs(b1))
          +0.5*obsleng[pn-1]*log(psig1/sigma[pn-1])-lambda*fabs(b2);
        rsq2=*rsq+0.5*obsleng[rn-1]*log(rsig2/sigma[rn-1])-lambda*fabs(b1)
          +0.5*obsleng[pn-1]*log(psig2/sigma[pn-1])+lambda*(fabs(cbeta2)-fabs(b2));

        if(rsq1<=rsq2)
        {
          *rsq=rsq1;
          if(dif1!=0)
          {
            beta[rindex]=cbeta1;
            norm[rn-1]=rnorm1;
            sigma[rn-1]=rsig1;
            for (k=(rn-1)*node;k<rn*node;k++)
            {
              inprod[k]-=fprod[(rn-1)*(node2-node)+(pn-1)*node+k]*dif1;
            }
          }
          if(b2!=0)
          {
            beta[pindex]=0.0;
            norm[pn-1]=pnorm1;
            sigma[pn-1]=psig1;
            dif2=-b2;
            for (k=(pn-1)*node;k<pn*node;k++)
            {
              inprod[k]-=fprod[(pn-1)*(node2-node)+(rn-1)*node+k]*dif2;
            }
          }
          if(cbeta1!=0)
          {
            G[rindex]=1;
            active[active_length]=i;
            active_length++;
          }
          else
            G[rindex]=0;
          G[pindex]=0;
          mad=Max(Max(fabs(dif1),fabs(b2)),mad);
        }
        else
        {
          *rsq=rsq2;
          if(b1!=0)
          {
            beta[rindex]=0.0;
            norm[rn-1]=rnorm2;
            sigma[rn-1]=rsig2;
            dif1=-b1;
            for (k=(rn-1)*node;k<rn*node;k++)
            {
              inprod[k]-=fprod[(rn-1)*(node2-node)+(pn-1)*node+k]*dif1;
            }
          }
          if(dif2!=0)
          {
            beta[pindex]=cbeta2;
            norm[pn-1]=pnorm2;
            sigma[pn-1]=psig2;
            for (k=(pn-1)*node;k<pn*node;k++)
            {
              inprod[k]-=fprod[(pn-1)*(node2-node)+(rn-1)*node+k]*dif2;
            }
          }
          if(cbeta2!=0)
          {
            G[pindex]=1;
            active[active_length]=i;
            active_length++;
          }
          else
            G[pindex]=0;
          G[rindex]=0;
          mad=Max(Max(fabs(b1),fabs(dif2)),mad);
        }
      }
    }

    if(mad<eps)
    {
      Free(active);
      break;
    }



    //active set convergence
    while(1)
    {
      mad=0.0;
      for(i=0;i<active_length;i++)
      {
        rn=eorder[active[i]*2];
        pn=eorder[active[i]*2+1];
        pindex=(pn-1)*node+rn-1;
        rindex=(rn-1)*node+pn-1;
        for (j=0;j<node2;j++)
          tG[j]=G[j];
        tG[rindex]=tG[pindex]=0;
        b1=beta[rindex];
        b2=beta[pindex];

        if(Cycle(node,tG,rn,pn))
        {
          denom1=fprod[(rn-1)*node2+(pn-1)*(node+1)];
          rnorm2=norm[rn-1]+2*b1*inprod[rindex]+denom1*b1*b1;
          absum1=inprod[rindex]+denom1*b1;
          tmpSigSq1=rnorm2/obsleng[rn-1];
          bda1=0.5*obsleng[rn-1]*denom1/(lambda*lambda);
          Delta1=bda1*bda1-(denom1*rnorm2-absum1*absum1)/(lambda*lambda);
          if(absum1>lambda*tmpSigSq1)
          {
            rsig1=bda1-sqrt(Delta1);
            cbeta1=(absum1-lambda*rsig1)/denom1;
          }
          else if(absum1<(-lambda*tmpSigSq1))
          {
            rsig1=bda1-sqrt(Delta1);
            cbeta1=(absum1+lambda*rsig1)/denom1;
          }
          else if(Delta1>=0 && fabs(absum1)>lambda*bda1)
          {
            root1=bda1+sqrt(Delta1);
            root2=bda1-sqrt(Delta1);
            root3=tmpSigSq1;
            if(absum1>0)
            {
              tmpCoef1=(absum1-lambda*root1)/denom1;
              tmpCoef2=(absum1-lambda*root2)/denom1;
              tmpCoef3=0;
            } else
            {
              tmpCoef1=(absum1+lambda*root1)/denom1;
              tmpCoef2=(absum1+lambda*root2)/denom1;
              tmpCoef3=0;
            }
            dLoss1=0.5*obsleng[rn-1]*log(root1)+lambda*fabs(tmpCoef1);
            dLoss2=0.5*obsleng[rn-1]*log(root2)+lambda*fabs(tmpCoef2);
            dLoss3=0.5*obsleng[rn-1]*log(root3)+lambda*fabs(tmpCoef3);
            rsig1= (dLoss1<dLoss2) ? ((dLoss1<dLoss3) ? root1 : root3) : ((dLoss2<dLoss3) ? root2 : root3);
            cbeta1= (dLoss1<dLoss2) ? ((dLoss1<dLoss3) ? tmpCoef1 : tmpCoef3) : ((dLoss2<dLoss3) ? tmpCoef2 : tmpCoef3);
          }
          else
          {
            rsig1=tmpSigSq1;
            cbeta1=0;
          }
          dif1=cbeta1-b1;
          if(dif1!=0)
          {
            rnorm1=norm[rn-1]-2*dif1*inprod[rindex]+denom1*dif1*dif1;
            *rsq=*rsq+0.5*obsleng[rn-1]*log(rsig1/sigma[rn-1])+lambda*(fabs(cbeta1)-fabs(b1));
            beta[rindex]=cbeta1;
            norm[rn-1]=rnorm1;
            sigma[rn-1]=rsig1;
            for (k=(rn-1)*node;k<rn*node;k++)
            {
              inprod[k]-=fprod[(rn-1)*(node2-node)+(pn-1)*node+k]*dif1;
            }
          }
          if(cbeta1!=0)
          {
            G[rindex]=1;
          }
          else
            G[rindex]=0;
          mad=Max(fabs(dif1),mad);
        }

        else if(Cycle(node,tG,pn,rn))
        {
          denom2=fprod[(pn-1)*node2+(rn-1)*(node+1)];
          pnorm1=norm[pn-1]+2*b2*inprod[pindex]+denom2*b2*b2;
          absum2=inprod[pindex]+denom2*b2;
          tmpSigSq2=pnorm1/obsleng[pn-1];
          bda2=0.5*obsleng[pn-1]*denom2/(lambda*lambda);
          Delta2=bda2*bda2-(denom2*pnorm1-absum2*absum2)/(lambda*lambda);
          if(absum2>lambda*tmpSigSq2)
          {
            psig2=bda2-sqrt(Delta2);
            cbeta2=(absum2-lambda*psig2)/denom2;
          }
          else if(absum2<(-lambda*tmpSigSq2))
          {
            psig2=bda2-sqrt(Delta2);
            cbeta2=(absum2+lambda*psig2)/denom2;
          }
          else if(Delta2>=0 && fabs(absum2)>lambda*bda2)
          {
            root1=bda2+sqrt(Delta2);
            root2=bda2-sqrt(Delta2);
            root3=tmpSigSq2;
            if(absum2>0)
            {
              tmpCoef1=(absum2-lambda*root1)/denom2;
              tmpCoef2=(absum2-lambda*root2)/denom2;
              tmpCoef3=0;
            } else
            {
              tmpCoef1=(absum2+lambda*root1)/denom2;
              tmpCoef2=(absum2+lambda*root2)/denom2;
              tmpCoef3=0;
            }
            dLoss1=0.5*obsleng[pn-1]*log(root1)+lambda*fabs(tmpCoef1);
            dLoss2=0.5*obsleng[pn-1]*log(root2)+lambda*fabs(tmpCoef2);
            dLoss3=0.5*obsleng[pn-1]*log(root3)+lambda*fabs(tmpCoef3);
            psig2= (dLoss1<dLoss2) ? ((dLoss1<dLoss3) ? root1 : root3) : ((dLoss2<dLoss3) ? root2 : root3);
            cbeta2= (dLoss1<dLoss2) ? ((dLoss1<dLoss3) ? tmpCoef1 : tmpCoef3) : ((dLoss2<dLoss3) ? tmpCoef2 : tmpCoef3);
          }
          else
          {
            psig2=tmpSigSq2;
            cbeta2=0;
          }
          dif2=cbeta2-b2;
          if(dif2!=0)
          {
            pnorm2=norm[pn-1]-2*dif2*inprod[pindex]+denom2*dif2*dif2;
            *rsq=*rsq+0.5*obsleng[pn-1]*log(psig2/sigma[pn-1])+lambda*(fabs(cbeta2)-fabs(b2));
            beta[pindex]=cbeta2;
            norm[pn-1]=pnorm2;
            sigma[pn-1]=psig2;
            for (k=(pn-1)*node;k<pn*node;k++)
            {
              inprod[k]-=fprod[(pn-1)*(node2-node)+(rn-1)*node+k]*dif2;
            }
          }
          if(cbeta2!=0)
          {
            G[pindex]=1;
          }
          else
            G[pindex]=0;
          mad=Max(fabs(dif2),mad);
        }

        else
        {
          denom1=fprod[(rn-1)*node2+(pn-1)*(node+1)];
          denom2=fprod[(pn-1)*node2+(rn-1)*(node+1)];

          rnorm2=norm[rn-1]+2*b1*inprod[rindex]+denom1*b1*b1;
          absum1=inprod[rindex]+denom1*b1;
          tmpSigSq1=rnorm2/obsleng[rn-1];
          bda1=0.5*obsleng[rn-1]*denom1/(lambda*lambda);
          Delta1=bda1*bda1-(denom1*rnorm2-absum1*absum1)/(lambda*lambda);
          if(absum1>lambda*tmpSigSq1)
          {
            rsig1=bda1-sqrt(Delta1);
            cbeta1=(absum1-lambda*rsig1)/denom1;
          }
          else if(absum1<(-lambda*tmpSigSq1))
          {
            rsig1=bda1-sqrt(Delta1);
            cbeta1=(absum1+lambda*rsig1)/denom1;
          }
          else if(Delta1>=0 && fabs(absum1)>lambda*bda1)
          {
            root1=bda1+sqrt(Delta1);
            root2=bda1-sqrt(Delta1);
            root3=tmpSigSq1;
            if(absum1>0)
            {
              tmpCoef1=(absum1-lambda*root1)/denom1;
              tmpCoef2=(absum1-lambda*root2)/denom1;
              tmpCoef3=0;
            } else
            {
              tmpCoef1=(absum1+lambda*root1)/denom1;
              tmpCoef2=(absum1+lambda*root2)/denom1;
              tmpCoef3=0;
            }
            dLoss1=0.5*obsleng[rn-1]*log(root1)+lambda*fabs(tmpCoef1);
            dLoss2=0.5*obsleng[rn-1]*log(root2)+lambda*fabs(tmpCoef2);
            dLoss3=0.5*obsleng[rn-1]*log(root3)+lambda*fabs(tmpCoef3);
            rsig1= (dLoss1<dLoss2) ? ((dLoss1<dLoss3) ? root1 : root3) : ((dLoss2<dLoss3) ? root2 : root3);
            cbeta1= (dLoss1<dLoss2) ? ((dLoss1<dLoss3) ? tmpCoef1 : tmpCoef3) : ((dLoss2<dLoss3) ? tmpCoef2 : tmpCoef3);
          }
          else
          {
            rsig1=tmpSigSq1;
            cbeta1=0;
          }
          dif1=cbeta1-b1;
          rnorm1=norm[rn-1]-2*dif1*inprod[rindex]+denom1*dif1*dif1;

          pnorm1=norm[pn-1]+2*b2*inprod[pindex]+denom2*b2*b2;
          absum2=inprod[pindex]+denom2*b2;
          tmpSigSq2=pnorm1/obsleng[pn-1];
          bda2=0.5*obsleng[pn-1]*denom2/(lambda*lambda);
          Delta2=bda2*bda2-(denom2*pnorm1-absum2*absum2)/(lambda*lambda);
          if(absum2>lambda*tmpSigSq2)
          {
            psig2=bda2-sqrt(Delta2);
            cbeta2=(absum2-lambda*psig2)/denom2;
          }
          else if(absum2<(-lambda*tmpSigSq2))
          {
            psig2=bda2-sqrt(Delta2);
            cbeta2=(absum2+lambda*psig2)/denom2;
          }
          else if(Delta2>=0 && fabs(absum2)>lambda*bda2)
          {
            root1=bda2+sqrt(Delta2);
            root2=bda2-sqrt(Delta2);
            root3=tmpSigSq2;
            if(absum2>0)
            {
              tmpCoef1=(absum2-lambda*root1)/denom2;
              tmpCoef2=(absum2-lambda*root2)/denom2;
              tmpCoef3=0;
            } else
            {
              tmpCoef1=(absum2+lambda*root1)/denom2;
              tmpCoef2=(absum2+lambda*root2)/denom2;
              tmpCoef3=0;
            }
            dLoss1=0.5*obsleng[pn-1]*log(root1)+lambda*fabs(tmpCoef1);
            dLoss2=0.5*obsleng[pn-1]*log(root2)+lambda*fabs(tmpCoef2);
            dLoss3=0.5*obsleng[pn-1]*log(root3)+lambda*fabs(tmpCoef3);
            psig2= (dLoss1<dLoss2) ? ((dLoss1<dLoss3) ? root1 : root3) : ((dLoss2<dLoss3) ? root2 : root3);
            cbeta2= (dLoss1<dLoss2) ? ((dLoss1<dLoss3) ? tmpCoef1 : tmpCoef3) : ((dLoss2<dLoss3) ? tmpCoef2 : tmpCoef3);
          }
          else
          {
            psig2=tmpSigSq2;
            cbeta2=0;
          }
          dif2=cbeta2-b2;
          pnorm2=norm[pn-1]-2*dif2*inprod[pindex]+denom2*dif2*dif2;

          psig1=pnorm1/obsleng[pn-1];
          rsig2=rnorm2/obsleng[rn-1];

          rsq1=*rsq+0.5*obsleng[rn-1]*log(rsig1/sigma[rn-1])+lambda*(fabs(cbeta1)-fabs(b1))
            +0.5*obsleng[pn-1]*log(psig1/sigma[pn-1])-lambda*fabs(b2);
          rsq2=*rsq+0.5*obsleng[rn-1]*log(rsig2/sigma[rn-1])-lambda*fabs(b1)
            +0.5*obsleng[pn-1]*log(psig2/sigma[pn-1])+lambda*(fabs(cbeta2)-fabs(b2));


          if(rsq1<=rsq2)
          {
            *rsq=rsq1;
            if(dif1!=0)
            {
              beta[rindex]=cbeta1;
              norm[rn-1]=rnorm1;
              sigma[rn-1]=rsig1;
              for (k=(rn-1)*node;k<rn*node;k++)
              {
                inprod[k]-=fprod[(rn-1)*(node2-node)+(pn-1)*node+k]*dif1;
              }
            }
            if(b2!=0)
            {
              beta[pindex]=0.0;
              norm[pn-1]=pnorm1;
              sigma[pn-1]=psig1;
              dif2=-b2;
              for (k=(pn-1)*node;k<pn*node;k++)
              {
                inprod[k]-=fprod[(pn-1)*(node2-node)+(rn-1)*node+k]*dif2;
              }
            }
            if(cbeta1!=0)
            {
              G[rindex]=1;
            }
            else
              G[rindex]=0;
            G[pindex]=0;
            mad=Max(Max(fabs(dif1),fabs(b2)),mad);
          }
          else
          {
            *rsq=rsq2;
            if(b1!=0)
            {
              beta[rindex]=0.0;
              norm[rn-1]=rnorm2;
              sigma[rn-1]=rsig2;
              dif1=-b1;
              for (k=(rn-1)*node;k<rn*node;k++)
              {
                inprod[k]-=fprod[(rn-1)*(node2-node)+(pn-1)*node+k]*dif1;
              }
            }
            if(dif2!=0)
            {
              beta[pindex]=cbeta2;
              norm[pn-1]=pnorm2;
              sigma[pn-1]=psig2;
              for (k=(pn-1)*node;k<pn*node;k++)
              {
                inprod[k]-=fprod[(pn-1)*(node2-node)+(rn-1)*node+k]*dif2;
              }
            }
            if(cbeta2!=0)
            {
              G[pindex]=1;
            }
            else
              G[pindex]=0;
            G[rindex]=0;
            mad=Max(Max(fabs(b1),fabs(dif2)),mad);
          }
        }
      }
      if(mad<eps)
        break;
    }
    Free(active);

  }
  Free(G);
  Free(tG);
  eor=NULL;
}



//coordinate descent on a grid of lambdas.
void CDGrid(int node,double *beta,double *lambda,double *inprod,double *fprod,int *eorder,double rsq,double eps,double *coef,double *RSQ,int *ITER,int nlam,double *sigma,double *sig,double *norm,int *obsleng,int *eorleng)
{
  int i,k,iter,node2=node*node;
  double bsum,*tmp=coef,*tmp2=sig;
  for (k=0;k<nlam+1;k++)
  {
    if(k>0)
    {
      bsum=0;
      for (i=0;i<node2;i++)
      {
        bsum+=fabs(beta[i]);
      }
      rsq+=(lambda[k]-lambda[k-1])*bsum;
    }
    iter=0;

    CoordinateDescent(node,beta,lambda[k],inprod,fprod,eorder,&rsq,eps,&iter,sigma,norm,obsleng,eorleng);

    for (i=0;i<node2;i++)
    {
      *tmp=beta[i];
      tmp++;
    }
    for (i=0;i<node;i++)
    {
      *tmp2=sigma[i];
      tmp2++;
    }
    RSQ[k]=rsq;
    ITER[k]=iter;
  }

  tmp=NULL;
  tmp2=NULL;
}

// a wrapper for CDGrid so the function can be called directly from R----Fei
void CDGridR(int *node,double *beta,double *lambda,double *inprod,double *fprod,int *eorder,double *rsq,double *eps,double *coef,double *RSQ,int *ITER,int *nlam,double *sigma,double *sig,double *norm,int *obsleng,int *eorleng)
{
  CDGrid(*node,beta,lambda,inprod,fprod,eorder,*rsq,*eps,coef,RSQ,ITER,*nlam,sigma,sig,norm,obsleng,eorleng);
}

// a wrapper for CDGrid so that we can export it using RcppExports----Jean
// [[Rcpp::export]]
SEXP gridCD(SEXP node,SEXP beta,SEXP lambda,SEXP inprod,SEXP fprod,SEXP eorder,SEXP rsq,SEXP eps,SEXP coef,SEXP RSQ,SEXP ITER,SEXP nlam,SEXP sigma,SEXP sig,SEXP norm,SEXP obsleng,SEXP eorleng)
{
  int xnode,*xeorder,*xITER,xnlam,*xobsleng,*xeorleng;
  double *xbeta,*xlambda,*xinprod,*xfprod,xrsq,xeps,*xcoef,*xRSQ,*xsigma,*xsig,*xnorm;

  xnode=INTEGER(node)[0];
  xeorder=INTEGER(eorder);
  xITER=INTEGER(ITER);
  xnlam=INTEGER(nlam)[0];
  xobsleng=INTEGER(obsleng);
  xeorleng=INTEGER(eorleng);
  xbeta=REAL(beta);
  xlambda=REAL(lambda);
  xinprod=REAL(inprod);
  xfprod=REAL(fprod);
  xrsq=REAL(rsq)[0];
  xeps=REAL(eps)[0];
  xcoef=REAL(coef);
  xRSQ=REAL(RSQ);
  xsigma=REAL(sigma);
  xsig=REAL(sig);
  xnorm=REAL(norm);

  CDGrid(xnode,xbeta,xlambda,xinprod,xfprod,xeorder,xrsq,xeps,xcoef,xRSQ,xITER,xnlam,xsigma,xsig,xnorm,xobsleng,xeorleng);

  return coef;
}




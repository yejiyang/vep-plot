/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * partdmesh.c
 *
 * This file reads in the element node connectivity array of a mesh and 
 * partitions both the elements and the nodes using KMETIS on the dual graph.
 *
 * Started 9/29/97
 * George
 *
 * $Id: partdmesh.c,v 1.1 1998/11/27 17:59:38 karypis Exp $
 *
 */

#include "metis.h"



/*************************************************************************
* Let the game begin
**************************************************************************/
//main(int argc, char *argv[])
void partdmain_(int *argc, char *argv1, char *argv2, int *arg, idxtype *epart1, idxtype *npart1,int *nee, int *npp)
{
  int i, j, ne, nn, etype, numflag=0, nparts, edgecut;
  idxtype *elmnts, *epart, *npart;

  timer IOTmr, DUALTmr;
  char etypestr[4][5] = {"TRI", "TET", "HEX", "QUAD"};
  GraphType graph;
  FILE *fp;

  /* argc=3; */
  if (*argc != 3) {
    printf("Usage: %s <meshfile> <nparts>\n",argv1);
    exit(0);
  }

  nparts = atoi(argv2);
  fp=fopen("part","r");
  if(fp==NULL)
  {
   printf("Open file %s Error.\n","part");
   exit (0);
  }
  fscanf(fp,"%d",&nparts);
  fclose(fp);
 /* printf ("argc===,%d %d\n",*argc,nparts);*/
  if (*arg != nparts) {
    printf("Fatal error, different partition information %d  %d \n",*arg,nparts); 
    exit(0);
  }
  
  if (nparts < 2) {
    printf("nparts must be greater than one.\n");
    exit(0);
  }
   
  cleartimer(IOTmr);
  cleartimer(DUALTmr);

  starttimer(IOTmr);
  elmnts = ReadMesh(argv1, &ne, &nn, &etype);
  stoptimer(IOTmr);

  epart = idxmalloc(ne, "main: epart");
  /* epart1 = idxmalloc(ne, "main: epart1");  */
  npart = idxmalloc(nn, "main: npart");
  /* npart1 = idxmalloc(nn, "main: npart1");  */
  *nee=ne;
  *npp=nn;

  printf("**********************************************************************\n");
//  printf("%s", METISTITLE);
//  printf("Mesh Information ----------------------------------------------------\n");
//  printf("  Name: %s, #Elements: %d, #Nodes: %d, Etype: %s\n\n", argv1, ne, nn, etypestr[etype-1]);
  printf("Partitioning Dual Graph... ------------------------------------------\n");
  starttimer(DUALTmr);
  METIS_PartMeshDual(&ne, &nn, elmnts, &etype, &numflag, &nparts, &edgecut, epart, npart);
  stoptimer(DUALTmr);
//  printf("epart[0]======%d\n",epart[0]);
//  printf("npart[0]======%d\n",npart[0]);
  for (i=0;i<ne;i++)
 {*(epart1+i)=epart[i];
 };
  for (i=0;i<nn;i++)
 {*(npart1+i)=npart[i];
 };
/*  printf ("asdasdadasdasdasdsasad %d   %d :\n",ne,ne); 
  for (i=0;i<ne;i++)
 	    {*(epart1+i)=epart[i];
			};
  for (i=0;i<ne;i++)
	    {*(npart1+i)=npart[i];
			};
  printf("aaaa%d %d %d %d\n",epart[100],epart1[100],npart[100],npart1[100]);
 
  printf("  %d-way Edge-Cut: %7d, Balance: %5.2f\n", nparts, edgecut, ComputeElementBalance(ne, nparts, epart));
*/
/*  starttimer(IOTmr);
  WriteMeshPartition(argv1, nparts, ne, epart, nn, npart);
  stoptimer(IOTmr);

*/
 /* printf("\nTiming Information --------------------------------------------------\n");
  printf("  I/O:          \t\t %7.3f\n", gettimer(IOTmr));
  printf("  Partitioning: \t\t %7.3f\n", gettimer(DUALTmr));
  printf("**********************************************************************\n");
*/
/*
  graph.nvtxs = nn;
  graph.xadj = idxmalloc(nn+1, "xadj");
  graph.vwgt = idxsmalloc(nn, 1, "vwgt");
  graph.adjncy = idxmalloc(20*nn, "adjncy");
  graph.adjwgt = idxsmalloc(20*nn, 1, "adjncy");

  METIS_MeshToNodal(&ne, &nn, elmnts, &etype, &numflag, graph.xadj, graph.adjncy);

  ComputePartitionInfo(&graph, nparts, npart);

  GKfree(&graph.xadj, &graph.adjncy, &graph.vwgt, &graph.adjwgt, LTERM);
*/

  GKfree(&elmnts, &epart, &npart, LTERM);

}



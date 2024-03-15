/******************************************************************************
 *	Function write_imgheader writes header information to a
 *	microstructure file
 *
 * 	Arguments:	file pointer
 * 				pointer to int xsize
 * 				pointer to int ysize
 * 				pointer to int zsize
 * 				pointer to float resolution
 *
 *	Returns:	int status flag (0 if okay, 1 if otherwise)
 *
 *	Programmer:	Jeffrey W. Bullard
 *				NIST
 *				100 Bureau Drive, Stop 8615
 *				Gaithersburg, Maryland  20899-8615
 *				USA
 *
 *				Phone:	301.975.5725
 *				Fax:	301.990.6891
 *				bullard@nist.gov
 *
 *	16 March 2004
 ******************************************************************************/
#include "auxmic.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int write_imgheader(FILE *fpout, int xsize, int ysize, int zsize, float res) {
  int status = 0;

  if (!fpout) {
    status = 1;
    return (status);
  }

  fprintf(fpout, "%s %s\n", VersionString, VersionNumber);
  fprintf(fpout, "%s %d\n", XSizeString, xsize);
  fprintf(fpout, "%s %d\n", YSizeString, ysize);
  fprintf(fpout, "%s %d\n", ZSizeString, zsize);
  fprintf(fpout, "%s %4.2f\n", ImgResString, res);

  return (status);
}

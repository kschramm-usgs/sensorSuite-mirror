/***************************************************************************
 * parse_one_file.c
 *
 * A program for libmseed parsing tests.
 *
 * Written by Chad Trabant, IRIS Data Management Center
 *
 * modified 2016.275
 ***************************************************************************/

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <libmseed.h>

#define PACKAGE "parse_one_file"
#define VERSION "[libmseed " LIBMSEED_VERSION " " PACKAGE " ]"

static flag verbose    = 0;
static flag ppackets   = 0;
static flag basicsum   = 0;
static flag tracegap   = 0;
static int printraw    = 0;
static int printdata   = 1;
static int reclen      = -1;
static char const * inputfile = "../test/TST5_00_BH0.512.seed";

static double timetol     = -1.0; /* Time tolerance for continuous traces */
static double sampratetol = -1.0; /* Sample rate tolerance for continuous traces */

int main (int argc, char **argv) {

  MSTraceList * mstl = NULL;
  MSRecord * msr     = NULL;

  int64_t totalrecs  = 0;
  int64_t totalsamps = 0;
  int retcode;

  /* Loop over the input file */
  while ((retcode = ms_readmsr (&msr, inputfile, reclen, NULL, NULL, 1,
                                printdata, verbose)) == MS_NOERROR)
  {
    totalrecs++;
    totalsamps += msr->samplecnt;

    if (tracegap)
    {
      mstl_addmsr (mstl, msr, 0, 1, timetol, sampratetol);
    }
    else
    {
      if ( printraw )
        ms_parse_raw (msr->record, msr->reclen, ppackets, -1);
      else
        msr_print (msr, ppackets);

      if (printdata && msr->numsamples > 0)
      {
        int line, col, cnt, samplesize;
        int lines = (msr->numsamples / 6) + 1;
        void *sptr;

        if ((samplesize = ms_samplesize (msr->sampletype)) == 0)
        {
          ms_log (2, "Unrecognized sample type: '%c'\n", msr->sampletype);
        }
        if (msr->sampletype == 'a')
        {
          char *ascii = (char *)msr->datasamples;
          int length  = msr->numsamples;

          ms_log (0, "ASCII Data:\n");

          /* Print maximum log message segments */
          while (length > (MAX_LOG_MSG_LENGTH - 1))
          {
            ms_log (0, "%.*s", (MAX_LOG_MSG_LENGTH - 1), ascii);
            ascii += MAX_LOG_MSG_LENGTH - 1;
            length -= MAX_LOG_MSG_LENGTH - 1;
          }

          /* Print any remaining ASCII and add a newline */
          if (length > 0)
          {
            ms_log (0, "%.*s\n", length, ascii);
          }
          else
          {
            ms_log (0, "\n");
          }
        }
        else
          for (cnt = 0, line = 0; line < lines; line++)
          {
            for (col = 0; col < 6; col++)
            {
              if (cnt < msr->numsamples)
              {
                sptr = (char *)msr->datasamples + (cnt * samplesize);

                if (msr->sampletype == 'i')
                  ms_log (0, "%10d  ", *(int32_t *)sptr);

                else if (msr->sampletype == 'f')
                  ms_log (0, "%10.8g  ", *(float *)sptr);

                else if (msr->sampletype == 'd')
                  ms_log (0, "%10.10g  ", *(double *)sptr);

                cnt++;
              }
            }
            ms_log (0, "\n");

            /* If only printing the first 6 samples break out here */
            if (printdata == 1)
              break;
          }
      }
    }
  }

  if (retcode != MS_ENDOFFILE)
    ms_log (2, "Cannot read %s: %s\n", inputfile, ms_errorstr (retcode));

  if (tracegap)
    mstl_printtracelist (mstl, 0, 1, 1);

  /* Make sure everything is cleaned up */
  ms_readmsr (&msr, NULL, 0, NULL, NULL, 0, 0, 0);

  if (mstl)
    mstl_free (&mstl, 0);

  if (basicsum)
    ms_log (1, "Records: %" PRId64 ", Samples: %" PRId64 "\n",
            totalrecs, totalsamps);

  return 0;
} /* End of main() */


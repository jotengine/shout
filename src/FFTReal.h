/*****************************************************************************
*                                                                            *
*       DIGITAL SIGNAL PROCESSING TOOLS                                      *
*       Version 1.03, 2001/06/15                                             *
*       (c) 1999 - Laurent de Soras                                          *
*                                                                            *
*       FFTReal.h                                                            *
*       Fourier transformation of real number arrays.                        *
*       Portable ISO C++                                                     *
*                                                                            *
* Tab = 3                                                                    *
*****************************************************************************/

#ifndef FFTReal_CURRENT_HEADER
#define FFTReal_CURRENT_HEADER


typedef float  flt_t;

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief The FFT procedures used by shout are written by Laurent de Soras
///
/// This is his legal message:
///
///  LEGAL
///
///  Source code may be freely used for any purpose, including commercial
///  applications. Programs must display in their "About" dialog-box (or
///  documentation) a text telling they use these routines by Laurent de Soras.
///  Modified source code can be distributed, but modifications must be clearly
///  indicated.
///
///  CONTACT
///
///  Laurent de Soras
///  92 avenue Albert 1er
///  92500 Rueil-Malmaison
///  France
///
///  ldesoras@club-internet.fr
/////////////////////////////////////////////////////////////////////////////////////////////////////

class FFTReal
{
  public:
    explicit  FFTReal (const long length);
             ~FFTReal ();
             
  public:
    void  do_fft (flt_t f [], const flt_t x []) const;
    void  do_ifft (const flt_t f [], flt_t x []) const;
    void  rescale (flt_t x []) const;

  private:
/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief The FFT procedures used by shout are written by Laurent de Soras
///
/// This is his legal message:
///
///  LEGAL
///
///  Source code may be freely used for any purpose, including commercial
///  applications. Programs must display in their "About" dialog-box (or
///  documentation) a text telling they use these routines by Laurent de Soras.
///  Modified source code can be distributed, but modifications must be clearly
///  indicated.
///
///  CONTACT
///
///  Laurent de Soras
///  92 avenue Albert 1er
///  92500 Rueil-Malmaison
///  France
///
///  ldesoras@club-internet.fr
/////////////////////////////////////////////////////////////////////////////////////////////////////
  /* Bit-reversed look-up table nested class */  
  class BitReversedLUT
  {
    public:
      explicit BitReversedLUT (const int nbr_bits);
              ~BitReversedLUT ();
      const long *  get_ptr () const
      {
        return (_ptr);
      }
    private:
      long * _ptr;
  };

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief The FFT procedures used by shout are written by Laurent de Soras
///
/// This is his legal message:
///
///  LEGAL
///
///  Source code may be freely used for any purpose, including commercial
///  applications. Programs must display in their "About" dialog-box (or
///  documentation) a text telling they use these routines by Laurent de Soras.
///  Modified source code can be distributed, but modifications must be clearly
///  indicated.
///
///  CONTACT
///
///  Laurent de Soras
///  92 avenue Albert 1er
///  92500 Rueil-Malmaison
///  France
///
///  ldesoras@club-internet.fr
/////////////////////////////////////////////////////////////////////////////////////////////////////
/* Trigonometric look-up table nested class */
  class  TrigoLUT
  {
    public:
      explicit      TrigoLUT (const int nbr_bits);
                   ~TrigoLUT ();
      const flt_t  *  get_ptr (const int level) const
      {
        return (_ptr + (1L << (level - 1)) - 4);
      };
    private:
      flt_t  *      _ptr;
  };

  const int             _nbr_bits;
  const long            _length;
  const BitReversedLUT  _bit_rev_lut;
  const TrigoLUT        _trigo_lut;
  const flt_t           _sqrt2_2;
  flt_t *               _buffer_ptr;


private:
  FFTReal (const FFTReal &other);
  const FFTReal&  operator = (const FFTReal &other);
  int            operator == (const FFTReal &other);
  int            operator != (const FFTReal &other);
};

#endif  // FFTReal_HEADER_INCLUDED


/*****************************************************************************

  LEGAL

  Source code may be freely used for any purpose, including commercial
  applications. Programs must display in their "About" dialog-box (or
  documentation) a text telling they use these routines by Laurent de Soras.
  Modified source code can be distributed, but modifications must be clearly
  indicated.

  CONTACT

  Laurent de Soras
  92 avenue Albert 1er
  92500 Rueil-Malmaison
  France

  ldesoras@club-internet.fr

*****************************************************************************/


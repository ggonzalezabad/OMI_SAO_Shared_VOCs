module hdfeos4_parameters
  !       Number type info codes

      integer   DFNT_FLOAT32, DFNT_FLOAT, DFNT_FLOAT64
      integer   DFNT_DOUBLE,  DFNT_FLOAT128

      parameter(DFNT_FLOAT32    = 5)
      parameter(DFNT_FLOAT      = 5)
      parameter(DFNT_FLOAT64    = 6)
      parameter(DFNT_DOUBLE     = 6)
      parameter(DFNT_FLOAT128   = 7)
 
      integer   DFNT_INT8,  DFNT_UINT8
      integer   DFNT_INT16, DFNT_UINT16
      integer   DFNT_INT32, DFNT_UINT32
      integer   DFNT_INT64, DFNT_UINT64
      integer   DFNT_INT128,DFNT_UINT128

      parameter(DFNT_INT8       = 20)
      parameter(DFNT_UINT8      = 21)
      parameter(DFNT_INT16      = 22)
      parameter(DFNT_UINT16     = 23)
      parameter(DFNT_INT32      = 24)
      parameter(DFNT_UINT32     = 25)
      parameter(DFNT_INT64      = 26)
      parameter(DFNT_UINT64     = 27)
      parameter(DFNT_INT128     = 28)
      parameter(DFNT_UINT128    = 29)

      integer   DFNT_UCHAR8, DFNT_UCHAR, DFNT_CHAR8
      integer   DFNT_CHAR,   DFNT_CHAR16, DFNT_UCHAR16

      parameter(DFNT_UCHAR8     = 3)
      parameter(DFNT_UCHAR      = 3)
      parameter(DFNT_CHAR8      = 4)
      parameter(DFNT_CHAR       = 4)
      parameter(DFNT_CHAR16     = 42)
      parameter(DFNT_UCHAR16    = 43)

      ! internal file access codes

      integer   DFACC_READ, DFACC_WRITE, DFACC_CREATE, DFACC_ALL
      integer   DFACC_RDONLY, DFACC_RDWR, DFACC_CLOBBER

      parameter(DFACC_READ               = 1)
      parameter(DFACC_WRITE              = 2)
      parameter(DFACC_CREATE             = 4)
      parameter(DFACC_ALL                = 7)
      parameter(DFACC_RDONLY             = 1)
      parameter(DFACC_RDWR               = 3)
      parameter(DFACC_CLOBBER            = 4)

      ! he4 functions

      integer  swdetach
      external swdetach
      integer  swclose
      external swclose
      
      integer   swopen 
      external  swopen
      integer   swattach
      external  swattach
      integer   swdiminfo
      external  swdiminfo
      integer   swrdattr
      external  swrdattr

      integer HDFE_NOMERGE
      parameter (HDFE_NOMERGE=0)
      integer HDFE_AUTOMERGE
      parameter (HDFE_AUTOMERGE=1)
      integer HDFE_MIDPOINT
      parameter (HDFE_MIDPOINT=0)
      integer HDFE_ENDPOINT
      parameter (HDFE_ENDPOINT=1)
      integer HDFE_INTERNAL
      parameter (HDFE_INTERNAL=0)
      integer HDFE_COMP_SKPHUFF
      parameter (HDFE_COMP_SKPHUFF=3)
      integer HDFE_COMP_NONE
      parameter (HDFE_COMP_NONE=0)
      integer HDFE_NOPREVSUB
      parameter (HDFE_NOPREVSUB=-1)

end module hdfeos4_parameters

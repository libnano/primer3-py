
#cdef extern from "libprimer3/thal.h":
#    ctypedef enum thal_alignment_type:
#        thal_any = 1
#        thal_end1 = 2
#        thal_end2 = 3
#        thal_hairpin = 4

#    ctypedef struct thal_args:
#        int debug # if non zero, print debugging info to stderr */
#        thal_alignment_type type
#        int maxLoop # maximum size of loop to consider; longer than 30 bp are not allowed */
#        double mv # concentration of monovalent cations */
#        double dv # concentration of divalent cations */
#        double dntp # concentration of dNTP-s */
#        double dna_conc # concentration of oligonucleotides */
#        double temp # temperature from which hairpin structures will be calculated */
#        int temponly # if non zero, print only temperature to stderr
#        int dimer   # if non zero, dimer structure is calculated
    
#    ctypedef struct thal_results:
#        char msg[255]
#        int no_structure # Added no structure (1 if no structure found)
#        double temp
#        double ds # Added entropy value
#        double dh # Added enthalpy value
#        double dg # Added gibbs free energy value
#        int align_end_1
#        int align_end_2



#    void thal(  const unsigned char*,
#                const unsigned char*,
#                const thal_args*,
#                thal_results*,
#                const int)

#    int get_thermodynamic_values(const char*, thal_results *)


cdef extern from "libprimer3/oligotm.h":
    ctypedef enum tm_method_type:
        breslauer_auto      = 0
        santalucia_auto     = 1

    ctypedef enum salt_correction_type:
        schildkraut    = 0
        santalucia     = 1
        owczarzy       = 2

    double seqtm(const char*,
                     double,
                     double,
                     double,
                     double,
                     int,
                     tm_method_type,
                     salt_correction_type,
                     )


from cpython.version cimport PY_MAJOR_VERSION

cdef unsigned char[:] _chars(s):
    cdef unsigned char[:] o
    if isinstance(s, str):
        # encode to the specific encoding used inside of the module
        o = memoryview(bytearray((<str>s).encode('utf8')))
        return o
    return memoryview(s)

cdef bytes _bytes(s):
    if PY_MAJOR_VERSION > 2 and isinstance(s, str):
        # encode to the specific encoding used inside of the module
        return  (<str>s).encode('utf8')
    else:
        return s

def loadThermoParams():
    cdef char*           param_path
    cdef thal_results    thalres
    import os
    PRIMER3_HOME = os.environ.get('PRIMER3HOME')
    ppath = os.path.join(PRIMER3_HOME, 'primer3_config/')
    ppathb = ppath.encode('utf-8')
    param_path = ppathb
    if get_thermodynamic_values(param_path, &thalres) != 0:
        raise IOError("couldn't load config file %s" % ppath)
loadThermoParams()

cdef class ThermoResult:
    #cdef thal_results thalres

    property temp:
        def __get__(self):
            return self.thalres.temp

    property ds:
        def __get__(self):
            return self.thalres.ds

    property dh:
        def __get__(self):
            return self.thalres.dh

    property dg:
        def __get__(self):
            return self.thalres.dg


    def __cinit__(self):
        self.thalres.no_structure = 0;
        self.thalres.ds = self.thalres.dh = self.thalres.dg = 0.0
        self.thalres.align_end_1 = self.thalres.align_end_2 = 0


cdef class ThermoAnalysis:
    #cdef thal_args thalargs
    
    # for melting temperature
    #cdef public int max_nn_length
    #cdef public int tm_method
    #cdef public int salt_correction_method

    def __cinit__(self, thal_type=1,
                    mv_conc=50,
                    dv_conc=0,
                    dntp_conc=0.8,
                    dna_conc=50,
                    temp_c=37,
                    max_loop=30,
                    temp_only=0,
                    debug=0,
                    max_nn_length=60,
                    tm_method=1,
                    salt_correction_method=1):
        self.thalargs.type = thal_type

        self.thalargs.mv = mv_conc;
        self.thalargs.dv = dv_conc;
        self.thalargs.dntp = dntp_conc;
        self.thalargs.dna_conc = dna_conc;
        self.thalargs.temp = temp_c + 273.15;
        self.thalargs.maxLoop = max_loop;
        self.thalargs.temponly = temp_only;
        self.thalargs.debug = debug;

        self.max_nn_length = max_nn_length;

        """
        'breslauer': 0,
        'santalucia': 1
        """
        self.tm_method = tm_method

        """
        'schildkraut': 0,
        'santalucia': 1,
        'owczarzy': 2
        """
        self.salt_correction_method = salt_correction_method
    # end def

    property thal_type:
        def __get__(self):
            return self.thalargs.type

    property mv_conc:
        def __get__(self):
            return self.thalargs.mv
        def __set__(self, value):
            self.thalargs.mv = value

    property dv_conc:
        def __get__(self):
            return self.thalargs.dv
        def __set__(self, value):
            self.thalargs.dv = value

    property dntp_conc:
        def __get__(self):
            return self.thalargs.dntp
        def __set__(self, value):
            self.thalargs.dntp_conc = value

    property dna_conc:
        def __get__(self):
            return self.thalargs.dna_conc
        def __set__(self, value):
            self.thalargs.dna_conc = value

    property temp:
        def __get__(self):
            return self.thalargs.temp - 273.15
        def __set__(self, value):
            self.thalargs.temp = value + 273.15

    property max_loop:
        def __get__(self):
            return self.thalargs.maxLoop

    cdef inline ThermoResult heterodimer_c(ThermoAnalysis self,
                                unsigned char*s1,
                                unsigned char* s2):
        cdef ThermoResult tr_obj = ThermoResult()

        self.thalargs.dimer = 1 
        self.thalargs.type = <thal_alignment_type> 1
        thal(<const unsigned char*> s1, <const unsigned char*> s2,
         <const thal_args *> &(self.thalargs), &(tr_obj.thalres), 0)
        return tr_obj
    # end def

    def heterodimer(ThermoAnalysis self, oligo1, oligo2):
        # first convert any unicode to a byte string and then
        # cooerce to a unsigned char *
        cdef unsigned char* s1 = <bytes> _bytes(oligo1)
        cdef unsigned char* s2 = <bytes> _bytes(oligo2)
        return ThermoAnalysis.heterodimer_c(<ThermoAnalysis> self, s1, s2)
    # end def

    cdef inline ThermoResult homodimer_c(ThermoAnalysis self, unsigned char*s1):
        cdef ThermoResult tr_obj = ThermoResult()

        self.thalargs.dimer = 1 
        self.thalargs.type = <thal_alignment_type> 1
        thal(<const unsigned char*> s1, <const unsigned char*> s1,
         <const thal_args *> &(self.thalargs), &(tr_obj.thalres), 0)
        return tr_obj
    # end def

    def homodimer(ThermoAnalysis self, oligo1):
        # first convert any unicode to a byte string and then
        # cooerce to a unsigned char *
        cdef unsigned char* s1 = <bytes> _bytes(oligo1)
        return ThermoAnalysis.homodimer_c(<ThermoAnalysis> self, s1)
    # end def

    cdef inline ThermoResult hairpin_c(ThermoAnalysis self, unsigned char*s1):
        cdef ThermoResult tr_obj = ThermoResult()

        self.thalargs.dimer = 1 
        self.thalargs.type = <thal_alignment_type> 4
        thal(<const unsigned char*> s1, <const unsigned char*> s1,
         <const thal_args *> &(self.thalargs), &(tr_obj.thalres), 0)
        return tr_obj
    # end def

    def hairpin(ThermoAnalysis self, oligo1):
        # first convert any unicode to a byte string and then
        # cooerce to a unsigned char *
        cdef unsigned char* s1 = <bytes> _bytes(oligo1)
        return ThermoAnalysis.hairpin_c(<ThermoAnalysis> self, s1)
    # end def


    cdef inline double meltingTemp_c(ThermoAnalysis self, char* s1):
        cdef thal_args *ta = &self.thalargs
        return seqtm(<const char*> s1, 
                            ta.dna_conc, 
                            ta.mv, 
                            ta.dv,
                            ta.dntp, 
                            self.max_nn_length, 
                            <tm_method_type>
                            self.tm_method,
                            <salt_correction_type> self.salt_correction_method)

    def meltingTemp(ThermoAnalysis self, oligo1):
        # first convert any unicode to a byte string and then
        # cooerce to a unsigned char *
        cdef char* s1 = <bytes> _bytes(oligo1)
        return ThermoAnalysis.meltingTemp_c(<ThermoAnalysis> self, s1)
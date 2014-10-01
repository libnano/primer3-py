
cdef extern from "libprimer3/thal.h":
    ctypedef enum thal_alignment_type:
        thal_any = 1
        thal_end1 = 2
        thal_end2 = 3
        thal_hairpin = 4

    ctypedef struct thal_args:
        int debug # if non zero, print debugging info to stderr */
        thal_alignment_type type
        int maxLoop # maximum size of loop to consider; longer than 30 bp are not allowed */
        double mv # concentration of monovalent cations */
        double dv # concentration of divalent cations */
        double dntp # concentration of dNTP-s */
        double dna_conc # concentration of oligonucleotides */
        double temp # temperature from which hairpin structures will be calculated */
        int temponly # if non zero, print only temperature to stderr
        int dimer   # if non zero, dimer structure is calculated
    
    ctypedef struct thal_results:
        char msg[255]
        int no_structure # Added no structure (1 if no structure found)
        double temp
        double ds # Added entropy value
        double dh # Added enthalpy value
        double dg # Added gibbs free energy value
        int align_end_1
        int align_end_2

    void thal(  const unsigned char *,
                const unsigned char *,
                const thal_args*,
                thal_results*,
                const int)

    int get_thermodynamic_values(const char*, thal_results *)

from cpython.version cimport PY_MAJOR_VERSION

cdef bytes _bytes(s):
    if PY_MAJOR_VERSION > 2 and isinstance(s, str):
        # encode to the specific encoding used inside of the module
        s = (<str>s).encode('utf8')
        return s
    else:
        return <bytes>s

def loadThermoParams():
    cdef char*           param_path
    cdef thal_results    thalres
    import os
    PRIMER3_HOME = os.environ.get('PRIMER3HOME')
    ppath = os.path.join(PRIMER3_HOME, 'primer3_config/')
    ppath = ppath.encode('utf-8')
    param_path = ppath
    if get_thermodynamic_values(param_path, &thalres) != 0:
        raise IOError("couldn't load config  file")
loadThermoParams()

cdef class ThermoResult:
    cdef thal_results thalres

    property temp:
        def __get__(self):
            return self.thalres.temp

    property dg:
        def __get__(self):
            return self.thalres.dg


    def __init__(self):
        self.thalres.no_structure = 0;
        self.thalres.ds = self.thalres.dh = self.thalres.dg = 0.0
        self.thalres.align_end_1 = self.thalres.align_end_2 = 0


cdef class ThermoAnalysis:
    cdef thal_args thalargs
    
    # for melting temperature
    cdef int max_nn_length
    cdef int tm_method
    cdef int salt_correction_method

    def __init__(self, thal_type=1,
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

    cdef inline heterodimer_c(ThermoAnalysis self, char* s1, char* s2):
        cdef ThermoResult tr_obj = ThermoResult()

        self.thalargs.dimer = 1 
        """
        this param is used by thal_main.c as a placeholder
        when determining whether the user has provided
        one or two sequences via the command line. if the
        user has only provided 1 sequence, then thal is
        run with the same sequence for oligo1 and oligo2
        """
        self.thalargs.type = <thal_alignment_type> 1
        thal(<const unsigned char*> s1, <const unsigned char*> s2,
         <const thal_args *> &(self.thalargs), &(tr_obj.thalres), 0)
        return tr_obj
    # end def

    def heterodimer(self, oligo1, oligo2):
        oligo1 = _bytes(oligo1)
        oligo2 = _bytes(oligo2)
        print("oligo1:", type(oligo1))
        print("oligo2:", type(oligo2))
        cdef char* s1 = oligo1
        cdef char* s2 = oligo2
        return ThermoAnalysis.heterodimer_c(<ThermoAnalysis> self, s1, s2)
    # end def

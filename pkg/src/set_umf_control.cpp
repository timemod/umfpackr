#include <Rcpp.h>
#include <umfpack.h>
#include <map>
using std::map;
using std::string;

static bool init = true;
static map<string, int> control_params;
static map<string, double> umfpack_strategy_opts;
static map<string, double> umfpack_ordering_opts;
static map<string, double> umfpack_scale_opts;

static void init_umf_control() {

    control_params["PRL"] = UMFPACK_PRL;
    control_params["DENSE_ROW"] = UMFPACK_DENSE_ROW;
    control_params["DENSE_COL"] = UMFPACK_DENSE_COL;
    control_params["PIVOT_TOLERANCE"] = UMFPACK_PIVOT_TOLERANCE;
    control_params["BLOCK_SIZE"] = UMFPACK_BLOCK_SIZE;
    control_params["STRATEGY`"] = UMFPACK_STRATEGY;
    control_params["ORDERING`"] = UMFPACK_ORDERING;
    control_params["ALLOC_INIT`"] = UMFPACK_ALLOC_INIT;
    control_params["IRSTEP`"] = UMFPACK_IRSTEP;
    control_params["FIXQ`"] = UMFPACK_FIXQ;
    control_params["AMD_DENSE"] = UMFPACK_AMD_DENSE;
    control_params["SYM_PIVOT_TOLERANCE"] = UMFPACK_SYM_PIVOT_TOLERANCE;
    control_params["SCALE"] = UMFPACK_SCALE;
    control_params["FRONT_ALLOC_INIT"] = UMFPACK_FRONT_ALLOC_INIT;
    control_params["DROPTOL"] = UMFPACK_DROPTOL;
    control_params["AGGRESSIVE"] = UMFPACK_AGGRESSIVE;
    control_params["SINGLETONS"] = UMFPACK_SINGLETONS;

    umfpack_strategy_opts["STRATEGY_AUTO"] = UMFPACK_STRATEGY_AUTO;
    umfpack_strategy_opts["STRATEGY_UNSYMMETRIC"] = UMFPACK_STRATEGY_UNSYMMETRIC;
    umfpack_strategy_opts["STRATEGY_SYMMETRIC"] = UMFPACK_STRATEGY_SYMMETRIC;

    umfpack_ordering_opts["ORDERING_CHOLMOD"] = UMFPACK_ORDERING_CHOLMOD;
    umfpack_ordering_opts["ORDERING_AMD"] = UMFPACK_ORDERING_AMD;
    umfpack_ordering_opts["ORDERING_GIVEN"] = UMFPACK_ORDERING_GIVEN;
    umfpack_ordering_opts["ORDERING_NONE"] = UMFPACK_ORDERING_NONE;
    umfpack_ordering_opts["ORDERING_METIS"] = UMFPACK_ORDERING_METIS;
    umfpack_ordering_opts["ORDERING_BEST"] = UMFPACK_ORDERING_BEST;
    umfpack_ordering_opts["ORDERING_USER"] = UMFPACK_ORDERING_USER;

    umfpack_scale_opts["SCALE_NONE"] = UMFPACK_SCALE_NONE;
    umfpack_scale_opts["SCALE_SUM"] = UMFPACK_SCALE_SUM;
    umfpack_scale_opts["SCALE_MAX"] = UMFPACK_SCALE_MAX;
}

double get_multiple_choice_option(const std::string option_name, 
                                 const Rcpp::CharacterVector option, 
                                 map<string, double> options) {
    std::string option_string = Rcpp::as<std::string>(option[0]);
    map<string, double>::iterator it;
    it = options.find(option_string);
    if (it != options.end()) {
        return it->second;
    } else {
       Rf_error("%s is not an UMFPACK %s option", option_string.c_str(), 
               option_name.c_str());
    }
}

void set_umf_control(double control[UMFPACK_CONTROL], Rcpp::List umf_control) {

    if (init) {
        init_umf_control();
        init = false;
    }

    // first set all umfpack default values
    umfpack_di_defaults(control);

    int nparam = umf_control.size();

    if (nparam == 0) return;

    Rcpp::CharacterVector param_names = umf_control.names();

    map<string, int>::iterator it;
    for (int i = 0; i < nparam; i++) {
        std::string name = Rcpp::as<std::string>(param_names[i]);
        if (name.compare("STRATEGY") == 0) {
            control[UMFPACK_STRATEGY] = 
                  get_multiple_choice_option(name, umf_control[i],
                                             umfpack_strategy_opts);
        } else if (name.compare("ORDERING") == 0) {
            control[UMFPACK_ORDERING] = 
                  get_multiple_choice_option(name, umf_control[i],
                                             umfpack_ordering_opts);
        } else if (name.compare("SCALE") == 0) {
            control[UMFPACK_SCALE] = 
                  get_multiple_choice_option(name, umf_control[i],
                                             umfpack_scale_opts);
        } else {
            it = control_params.find(name);
            if (it != control_params.end()) {
                int index = it->second;
                SEXP option = umf_control[i];
                switch( TYPEOF(option)) {
                case REALSXP:
                case INTSXP:
                case LGLSXP:
                    {
                    Rcpp::NumericVector option = umf_control[i];
                    control[index] = option[0];
                    break;
                    }
                default:
                    Rf_error("%s is not an numeric or logical value", name.c_str());
                }
            } else {
                Rf_error("%s is not an UMFPACK parameter", name.c_str());
            }
       }
    }
}

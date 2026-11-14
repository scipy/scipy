#ifndef IPX_CONTROL_H_
#define IPX_CONTROL_H_

#include <fstream>
#include <ostream>
#include <sstream>
#include <string>
#include "io/HighsIO.h"
#include "ipm/ipx/ipx_internal.h"
#include "ipm/ipx/multistream.h"
#include "ipm/ipx/timer.h"

namespace ipx {

// Class Control handles
// (1) accessing user parameters,
// (2) solver output,
// (3) solver interruption.
// Currently the solver is interrupted only by time limit. The interrupt
// mechanism could also be used if we run IPX and a simplex code concurrently
// and want to interrupt IPX when the simplex finished. For that reason a
// Control object cannot be copied; assuming that one thread sets an interrupt
// signal (not implemented yet), a call to control.InterruptCheck() from any
// part of the solver must return nonzero. Hence we must only have references or
// pointers to a single Control object in the whole of IPX.

class Control {
public:
    Control();

    // disallow copy and copy construction
    Control& operator=(const Control&) = delete;
    Control(const Control&) = delete;

    // disallow move and move construction
    Control& operator=(Control&&) = delete;
    Control(const Control&&) = delete;

    // Returns IPX_ERROR_* if interrupt is requested, 0 otherwise.
    Int InterruptCheck(const Int ipm_iteration_count = -1) const;

    // Returns output streams for log and debugging messages. The streams
    // evaluate to false if they discard output, so that we can write
    //
    //   if (control.Debug(3))
    //     control.Debug(3) << expensive_computation(...) << '\n';
    //
    // If the debug level is < 3, expensive_computation() is not performed.
    void hLog(std::stringstream& logging) const;
    void hLog(const std::string& str) const;
    void hDebug(std::stringstream& logging, Int level=1) const;
    std::ostream& Debug(Int level=1) const;

    // Sends logging to HiGHS logging or the log stream according to
    // parameters.highs_logging, if >= parameters.print_interval
    // seconds have been elapsed since the last call to IntervalLog()
    // or to ResetPrintInterval().
    void hIntervalLog(std::stringstream& logging) const;
    void ResetPrintInterval() const;

    double Elapsed() const;     // total runtime

    ipxint dualize() const { return parameters_.dualize; }
    ipxint scale() const { return parameters_.scale; }
    ipxint ipm_maxiter() const { return parameters_.ipm_maxiter; }
    double ipm_feasibility_tol() const {
        return parameters_.ipm_feasibility_tol; }
    double ipm_optimality_tol() const { return parameters_.ipm_optimality_tol; }
    double ipm_drop_primal() const { return parameters_.ipm_drop_primal; }
    double ipm_drop_dual() const { return parameters_.ipm_drop_dual; }
    double kkt_tol() const { return parameters_.kkt_tol; }
    ipxint crash_basis() const { return parameters_.crash_basis; }
    double dependency_tol() const { return parameters_.dependency_tol; }
    double volume_tol() const { return parameters_.volume_tol; }
    ipxint rows_per_slice() const { return parameters_.rows_per_slice; }
    ipxint maxskip_updates() const { return parameters_.maxskip_updates; }
    ipxint lu_kernel() const { return parameters_.lu_kernel; }
    double lu_pivottol() const { return parameters_.lu_pivottol; }
    ipxint run_crossover() const { return parameters_.run_crossover; }
    double start_crossover_tol() const { return parameters_.start_crossover_tol; }
    double pfeasibility_tol() const { return parameters_.pfeasibility_tol; }
    double dfeasibility_tol() const { return parameters_.dfeasibility_tol; }
    ipxint switchiter() const { return parameters_.switchiter; }
    ipxint stop_at_switch() const { return parameters_.stop_at_switch; }
    ipxint update_heuristic() const { return parameters_.update_heuristic; }
    ipxint maxpasses() const { return parameters_.maxpasses; }
    bool reportBasisData() const { return parameters_.analyse_basis_data; }
    ipxint runCentring() const{return parameters_.run_centring; }
    ipxint maxCentringSteps() const{return parameters_.max_centring_steps; }
    double centringRatioTolerance() const{return parameters_.centring_ratio_tolerance; }
    double centringRatioReduction() const {return parameters_.centring_ratio_reduction; }
    double centringAlphaScaling() const{return parameters_.centring_alpha_scaling; }
    ipxint badProductsTolerance() const{return parameters_.bad_products_tolerance; }
    bool timelessLog() const{return parameters_.timeless_log; }
    double timeLimit() const {return parameters_.time_limit; }
    void setTimerOffset(const double offset) { this->timer_.offset_ = offset; }


    const Parameters& parameters() const;
    void parameters(const Parameters& new_parameters);
    void callback(HighsCallback* callback);

    // Opens the log file defined in parameters.logfile, if any.
    // Ignores if an error occurs; in this case no file log is written.
    void OpenLogfile();

    // Closes an opened log file, if any.
    void CloseLogfile();

    // Resets the total runtime counter.
    void ResetTimer();

private:
    void MakeStream();           // composes output_
    Parameters parameters_;
    HighsCallback* callback_ = nullptr;
    std::ofstream logfile_;
    Timer timer_;                // total runtime
    mutable Timer interval_;     // time since last interval log
    mutable Multistream output_; // forwards to logfile and/or console
    mutable Multistream dummy_;  // discards everything
};

// Formats integer, string literal or floating point value into a string of
// size at least @width. So we can put formatted output into a stream without
// altering the stream's state.
std::string Format(Int i, int width);
std::string Format(const char* c, int width);
std::string Format(double d, int width, int prec,
                   std::ios_base::fmtflags floatfield);

// shortcuts
inline std::string Scientific(double d, int width, int prec) {
    return Format(d, width, prec, std::ios_base::scientific);
}
inline std::string Fixed(double d, int width, int prec) {
    return Format(d, width, prec, std::ios_base::fixed);
}
inline std::string sci2(double d) { return Scientific(d,0,2); }
inline std::string sci8(double d) { return Scientific(d,0,8); }
inline std::string fix2(double d) { return Fixed(d,0,2); }
inline std::string fix8(double d) { return Fixed(d,0,8); }
inline std::string time(double d) { return Fixed(d,8,1); }

// Formats @text into a line of fixed width and indentation.
// This is used to print messages like
//
//     Number of variables:                                1464
//     Number of constraints:                              696
//
// consistently via control.hLog(h_logging_stream) using 
//
//   h_logging_stream << Textline("Number of variables:") << 1464 << '\n'
//                    << Textline("Number of constraints:") << 696 << '\n';
//
template <typename T>
std::string Textline(const T& text) {
    std::ostringstream s;
    s << "    ";
    s.setf(std::ios_base::left);
    s.width(52);
    s << text;
    return s.str();
}

}  // namespace ipx

#endif  // IPX_CONTROL_H_

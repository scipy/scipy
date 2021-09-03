/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsTimer.h
 * @brief Profiling facility for computational components in HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef UTIL_HIGHSTIMER_H_
#define UTIL_HIGHSTIMER_H_

#include <cassert>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>

/**
 * @brief Clock record structure
 */
/*
struct HighsClockRecord {
  int calls;
  double start;
  double time;
  std::string name;
  std::string ch3_name;
};
*/
/**
 * @brief Class for profiling facility for computational components in HiGHS
 */
class HighsTimer {
 public:
  HighsTimer() {
    num_clock = 0;
    int i_clock = clock_def("Run HiGHS", "RnH");
    assert(i_clock == 0);
    run_highs_clock = i_clock;

    presolve_clock = clock_def("Presolve", "Pre");
    solve_clock = clock_def("Solve", "Slv");
    postsolve_clock = clock_def("Postsolve", "Pst");
  }

  /**
   * @brief Define a clock
   */
  int clock_def(
      const char* name,  //!< Full-length name (<=16 characters) for the clock
      const char* ch3_name  //!< 3-character name for the clock
  ) {
    int i_clock = num_clock;
    clock_num_call.push_back(0);
    clock_start.push_back(initial_clock_start);
    clock_time.push_back(0);
    clock_names.push_back(name);
    clock_ch3_names.push_back(ch3_name);
    num_clock++;
    return i_clock;
  }

  /**
   * @brief Zero an external clock record
   */
  /*
    void clockInit(HighsClockRecord& x_clock  //!< Record for the external clock
    ) {
      x_clock.calls = 0;
      x_clock.start = 0;
      x_clock.time = 0;
      x_clock.name = "";
      x_clock.ch3_name = "";
    }
    */

  /**
   * @brief Add to an external clock record
   */
  /*
    void clockAdd(HighsClockRecord x_clock,  //!< Record for the external clock
                  int i_clock                //!< Clock of record to be added
    ) {
      assert(i_clock >= 0);
      assert(i_clock < num_clock);
      x_clock.calls += clock_num_call[i_clock];
      x_clock.start = initial_clock_start;
      x_clock.time += clock_time[i_clock];
    }
    */

  /**
   * @brief Reset a HighsTimer instance to its state after the
   * constructor
   */
  void resetHighsTimer() {
    this->num_clock = 0;
    this->clock_num_call.clear();
    this->clock_start.clear();
    this->clock_names.clear();
    this->clock_ch3_names.clear();
    int i_clock = clock_def("Run HiGHS", "RnH");
    assert(i_clock == 0);
    this->run_highs_clock = i_clock;
    this->presolve_clock = clock_def("Presolve", "Pre");
    this->solve_clock = clock_def("Solve", "Slv");
    this->postsolve_clock = clock_def("Postsolve", "Pst");
  }

  /**
   * @brief Zero the data for all clocks
   */
  void zeroAllClocks() {
    for (int i = 0; i < num_clock; i++) {
      clock_num_call[i] = 0;
      clock_start[i] = initial_clock_start;
      clock_time[i] = 0;
    }
  }

  /**
   * @brief Start a clock
   */
  void start(int i_clock  //!< Index of the clock to be started
  ) {
    assert(i_clock >= 0);
    assert(i_clock < num_clock);
    // Check that the clock's been stopped. It should be set to
    // getWallTime() >= 0 (or initialised to initial_clock_start > 0)
#ifdef HiGHSDEV
    if (clock_start[i_clock] <= 0) {
      printf("recordStart [%2d] (%s) is %11.4g: _num_call = %d\n", i_clock,
             clock_names[i_clock].c_str(), clock_start[i_clock],
             clock_num_call[i_clock]);
      fflush(stdout);
    }
#endif
    assert(clock_start[i_clock] > 0);
    // Set the start to be the negation of the WallTick to check that
    // the clock's been started when it's next stopped
    clock_start[i_clock] = -getWallTime();
  }

  /**
   * @brief Stop a clock
   */
  void stop(int i_clock  //!< Index of the clock to be stopped
  ) {
    assert(i_clock >= 0);
    assert(i_clock < num_clock);
    // Check that the clock's been started. It should be set to
    // -getWallTime() <= 0
#ifdef HiGHSDEV
    if (clock_start[i_clock] > 0) {
      printf("recordFinish[%2d] (%s) is %11.4g: _num_call = %d\n", i_clock,
             clock_names[i_clock].c_str(), clock_start[i_clock],
             clock_num_call[i_clock]);
      fflush(stdout);
    }
#endif
    assert(clock_start[i_clock] < 0);
    double wall_time = getWallTime();
    double callClockTimes = wall_time + clock_start[i_clock];
    clock_time[i_clock] += callClockTimes;
    clock_num_call[i_clock]++;
    // Set the start to be the WallTick to check that the clock's been
    // stopped when it's next started
    clock_start[i_clock] = wall_time;
  }

  /**
   * @brief Read the time of a clock
   */
  double read(int i_clock  //!< Index of the clock to be read
  ) {
    assert(i_clock >= 0);
    assert(i_clock < num_clock);
    double read_time;
    if (clock_start[i_clock] < 0) {
      // The clock's been started, so find the current time
      double wall_time = getWallTime();
      read_time = clock_time[i_clock] + wall_time + clock_start[i_clock];
    } else {
      // The clock is currently stopped, so read the current time
      read_time = clock_time[i_clock];
    }
    return read_time;
  }

  /**
   * @brief Start the RunHighs clock
   */
  void startRunHighsClock() { start(run_highs_clock); }

  /**
   * @brief Stop the RunHighs clock
   */
  void stopRunHighsClock() { stop(run_highs_clock); }

  /**
   * @brief Read the RunHighs clock
   */
  double readRunHighsClock() { return read(run_highs_clock); }

  /**
   * @brief Test whether the RunHighs clock is running
   */
  bool runningRunHighsClock() { return clock_start[run_highs_clock] < 0; }

  /**
   * @brief Report timing information for the clock indices in the list
   */
  void report(const char* grep_stamp,  //!< Character string used to extract
                                       //!< output using grep
              std::vector<int>& clock_list,  //!< List of indices to report
              double ideal_sum_time = 0  //!< Ideal value for times to sum to
  ) {
    double tl_per_cent_report = 1.0;
    report_tl(grep_stamp, clock_list, ideal_sum_time, tl_per_cent_report);
  }

  void report_tl(
      const char*
          grep_stamp,  //!< Character string used to extract output using grep
      std::vector<int>& clock_list,  //!< List of indices to report
      double ideal_sum_time = 0,     //!< Ideal value for times to sum to
      double tl_per_cent_report =
          0  //!< Lower bound on percentage of total time
             //!< before an individual clock is reported
  ) {
    int num_clock_list_entries = clock_list.size();

    // Check validity of the clock list and check no clocks are still running
    for (int i = 0; i < num_clock_list_entries; i++) {
      assert(clock_list[i] >= 0);
      assert(clock_list[i] < num_clock);
      // Check that the clock's not still running. It should be set to
      // getWallTime() >= 0 (or initialised to initial_clock_start > 0)
#ifdef HiGHSDEV
      int i_clock = clock_list[i];
      if (clock_start[i_clock] <= 0) {
        printf(
            "Clock %2d (%s) is still running: Start = %11.4g: "
            "_num_call = %d\n",
            i_clock, clock_names[i_clock].c_str(), clock_start[i_clock],
            clock_num_call[i_clock]);
        fflush(stdout);
      }
#endif
      assert(clock_start[clock_list[i]] > 0);
    }
    // Determine whether there are any times to report
    int sum_calls = 0;
    for (int i = 0; i < num_clock_list_entries; i++)
      sum_calls += clock_num_call[clock_list[i]];
    if (!sum_calls) return;

    // Report in one line the per-mille contribution from each clock
    // First give the 3-character clock names as column headers
    printf("%s-name  ", grep_stamp);
    for (int i = 0; i < num_clock_list_entries; i++) {
      printf(" %-3s", clock_ch3_names[clock_list[i]].c_str());
    }
    printf("\n");

    // Then give the per-mille contribution relative to the total
    // HiGHS run time, and then relative to the sum of times for these
    // clocks
    double current_run_highs_time = readRunHighsClock();
    double sum_clock_times = 0;
    for (int passNum = 0; passNum < 3; passNum++) {
      // Don't write out if there's no ideal time
      if (passNum == 1 && ideal_sum_time <= 0) continue;
      double suPerMille = 0;
      if (passNum == 0) {
        printf("%s-total ", grep_stamp);
      } else if (passNum == 1) {
        printf("%s-ideal ", grep_stamp);
      } else {
        printf("%s-local ", grep_stamp);
      }
      for (int i = 0; i < num_clock_list_entries; i++) {
        int i_clock = clock_list[i];
        double perMille;
        if (passNum == 0) {
          perMille = 1000.0 * clock_time[i_clock] / current_run_highs_time;
        } else if (passNum == 1) {
          perMille = 1000.0 * clock_time[i_clock] / ideal_sum_time;
        } else {
          perMille = 1000.0 * clock_time[i_clock] / sum_clock_times;
        }
        int int_PerMille = (perMille + 0.5);  // Forcing proper rounding
        if (int_PerMille > 0) {
          printf("%4d", int_PerMille);  // Just in case one time is 1000!
        } else {
          printf("    ");  // Just in case one time is 1000!
        }
        suPerMille += perMille;
        if (passNum == 0) {
          sum_clock_times += clock_time[i_clock];
        }
      }
      int int_sum_permille = (suPerMille + 0.5);  // Forcing proper rounding
      printf(" per mille: Sum = %4d", int_sum_permille);
      printf("\n");
    }

    // Report one line per clock, the time, number of calls and time per call
    printf("%s-time  Operation         :    Time     ( Total", grep_stamp);
    if (ideal_sum_time > 0) printf(";  Ideal");
    printf(";  Local):    Calls  Time/Call\n");
    // Convert approximate seconds
    double sum_time = 0;
    for (int i = 0; i < num_clock_list_entries; i++) {
      int i_clock = clock_list[i];
      double time = clock_time[i_clock];
      double percent_run_highs = 100.0 * time / current_run_highs_time;
      double percent_sum_clock_times = 100.0 * time / sum_clock_times;
      double time_per_call = 0;
      if (clock_num_call[i_clock] > 0) {
        time_per_call = time / clock_num_call[i_clock];
        if (percent_sum_clock_times >= tl_per_cent_report) {
          printf("%s-time  %-18s: %11.4e (%5.1f%%", grep_stamp,
                 clock_names[i_clock].c_str(), time, percent_run_highs);
          if (ideal_sum_time > 0) {
            double percent_ideal = 100.0 * time / ideal_sum_time;
            printf("; %5.1f%%", percent_ideal);
          }
          printf("; %5.1f%%):%9d %11.4e\n", percent_sum_clock_times,
                 clock_num_call[i_clock], time_per_call);
        }
      }
      sum_time += time;
    }
    double percent_sum_clock_times = 100.0;
    double percent_run_highs = 100.0 * sum_time / current_run_highs_time;
    printf("%s-time  SUM               : %11.4e (%5.1f%%", grep_stamp, sum_time,
           percent_run_highs);
    if (ideal_sum_time > 0) {
      double percent_ideal = 100.0 * sum_time / ideal_sum_time;
      printf("; %5.1f%%", percent_ideal);
    }
    printf("; %5.1f%%)\n", percent_sum_clock_times);
    printf("%s-time  TOTAL             : %11.4e\n", grep_stamp,
           current_run_highs_time);
  }

  /**
   * @brief Return the current wall-clock time
   */
  double getWallTime() {
    using namespace std::chrono;
    return duration_cast<duration<double> >(
               wall_clock::now().time_since_epoch())
        .count();
  }

  virtual ~HighsTimer() = default;

  // private:
  using wall_clock = std::chrono::high_resolution_clock;
  using time_point = wall_clock::time_point;

  // Dummy positive start time for clocks - so they can be checked as
  // having been stopped
  const double initial_clock_start = 1.0;

  int num_clock = 0;
  std::vector<int> clock_num_call;
  std::vector<double> clock_start;
  std::vector<double> clock_time;
  std::vector<std::string> clock_names;
  std::vector<std::string> clock_ch3_names;
  // The index of the RunHighsClock - should always be 0
  int run_highs_clock;
  // Fundamental Highs clocks
  int presolve_clock;
  int solve_clock;
  int postsolve_clock;
};

#endif /* UTIL_HIGHSTIMER_H_ */

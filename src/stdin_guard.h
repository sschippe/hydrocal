// SPDX-License-Identifier: MIT
/**
 * @file stdin_guard.h
 * @brief End of input ends hydrocal instead of making it loop forever.
 *
 * hydrocal reads all its input from stdin, with scanf and with cin, at
 * about 400 places.  Almost none of them check for the end of the input,
 * so a script that stops in the middle of a prompt leaves the variables
 * unchanged, and every loop that asks again then runs at full speed and
 * writes prompts without end.
 *
 * Including this header in a source file that reads input does two things:
 *  - scanf(...) becomes a call of hydrocal_scanf, which behaves like scanf
 *    but quits the program when it returns EOF;
 *  - hydrocal_install_stdin_guard() puts a stream buffer on cin that quits
 *    the program when cin reaches the end of the input.
 *
 * Text that is not a number is not end of input and is handled as before.
 * fscanf and sscanf are not affected.
 */
#ifndef HYDROCAL_STDIN_GUARD_H
#define HYDROCAL_STDIN_GUARD_H

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <streambuf>

inline void hydrocal_quit_at_end_of_input(void) {
  std::fflush(stdout);
  std::printf("\n\n End of input, hydrocal quits.\n");
  std::fflush(stdout);
  std::exit(0);
}

/// scanf that quits the program at the end of the input
inline int hydrocal_scanf(const char *format, ...) {
  va_list args;
  va_start(args, format);
  int nread = std::vscanf(format, args);
  va_end(args);
  if (nread == EOF)
    hydrocal_quit_at_end_of_input();
  return nread;
}

/// stream buffer for cin that quits the program at the end of the input
class HydrocalEofGuard : public std::streambuf {
public:
  explicit HydrocalEofGuard(std::streambuf *source) : source_(source) {}

protected:
  int_type underflow() override {
    int_type c = source_->sgetc();
    if (traits_type::eq_int_type(c, traits_type::eof()))
      hydrocal_quit_at_end_of_input();
    return c;
  }
  int_type uflow() override {
    int_type c = source_->sbumpc();
    if (traits_type::eq_int_type(c, traits_type::eof()))
      hydrocal_quit_at_end_of_input();
    return c;
  }

private:
  std::streambuf *source_;
};

/// call once at the start of main, before any input is read
inline void hydrocal_install_stdin_guard(void) {
  static HydrocalEofGuard guard(std::cin.rdbuf());
  std::cin.rdbuf(&guard);
}

#ifndef HYDROCAL_NO_SCANF_GUARD
#define scanf hydrocal_scanf
#endif

#endif // HYDROCAL_STDIN_GUARD_H

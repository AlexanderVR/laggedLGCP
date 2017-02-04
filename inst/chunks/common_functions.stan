vector interleave(int N, vector x0, vector x1);
vector fft_sym(int N, vector x, vector y);
vector ifft_shift(int N, vector x, vector y, real h);
vector ifft_shift_centered(int M, vector x, vector y, real h);
vector matern_eigen(int M_nonzero, int N, real nu, real lengthscale, real delta);
vector transform_to_matern(vector re, vector im, int N, real nu, real lengthscale, real bin_width, real shift);
vector ifft_shift(int N, vector x, vector y, real h) {
  // compute inverse FFT of a complex conjugate-symmetric signal.
  // Also, apply an interpolatory shift by h * delta,
  // where delta is the sampling period.
  // Returns real vector of length N.

  vector[N / 2 + 1] x_new;
  vector[N / 2 - 1] y_new;
  x_new = x;

  // Do the shift by h bins.
  // For Inverse FFT conjugate the input signal
  for (j in 1:(N / 2 - 1)) {
    x_new[j + 1] = x[j + 1] * cos(2*pi() * j * h / N) - y[j] * sin(2*pi() * j * h /N);
    y_new[j] = -y[j] * cos(2*pi() * j * h / N) - x[j + 1] * sin(2*pi() * j * h / N);
  }
  return fft_sym(N, x_new, y_new);
}

vector ifft_shift_centered(int M, vector x, vector y, real h) {
  // compute inverse FFT of a complex conjugate-symmetric signal X
  // with zero DC and Nyquist components (X_0 = X_{N/2} = 0).
  // Also, apply an interpolatory shift by h * delta,
  // where delta is the sampling period.
  // Returns real vector of length N.

  vector[M + 2] x_new;
  vector[M] y_new;
  int N;

  N = 2 * (M + 1);
  x_new = rep_vector(0, M + 2);
  y_new = rep_vector(0, M);

  // Do the shift by h bins.
  // For Inverse FFT conjugate the input signal
  for (j in 1:num_elements(x)) {
    x_new[j + 1] = x[j] * cos(2*pi() * j * h / N) - y[j] * sin(2*pi() * j * h / N);
    y_new[j] = -y[j] * cos(2*pi() * j * h / N) - x[j] * sin(2*pi() * j * h / N);
  }
  return fft_sym(N, x_new, y_new);
}


vector interleave(int N, vector x0, vector x1) {
  // interlace elements of x0 and x1 into a vector of size 2N
  vector[2 * N] x;
  for (n in 1:N) {
    x[2 * n - 1] = x0[n];
    x[2 * n]     = x1[n];
  }
  return x;
}

vector fft_sym(int N, vector x, vector y) {
  // compute Fast Fourier Transform of a complex conjugate-symmetric signal
  // Conj(X[k]) = X[N - k], where N = 2^(j+1). Such signals have real FFT.
  // X is thus described by Real(X[0 : 2^j]) and Imag(X[1 : 2^j - 1])
  // Inputs: x must be of size N/2 + 1, y of size N/2 - 1.
  // Output: a real vector of length N.
  // Uses a "decimation-in-frequency" algorithm, which allows all
  // computations to stay real without wasting effort on computing the
  // (0-valued) complex entries of FFT(X)

  // The main recursion:
  // fft_sym(x, y) |--> interleave(fft_sym(x0, y0), fft_sym(x1, y1))
  // where xi, yi are of size N/4 + 1 and N/4 - 1, respectively.

  int M;
  // the real coefficients resulting from one step of the alg.
  vector[N / 4 + 1] x0;
  vector[N / 4 + 1] x1;

  M = N / 4;
  // Do one step of the decimation-in-frequency algorithm
  for (m in 0 : M) {
    x0[m + 1] =  x[m + 1] + x[2*M - m + 1];
    x1[m + 1] = (x[m + 1] - x[2*M - m + 1]) * cos(2*pi() * m / N);
  }

  for (m in 1 : M) {
    x1[m + 1] = x1[m + 1] + (y[m] + y[2*M - m]) * sin(2*pi() * m / N);
  }

  if (N == 4) {
    // Don't bother to compute the complex coefficients -- they are zero!
    vector[4] xout;
    xout[1] = x0[1] + x0[2]; xout[2] = x1[1] + x1[2]; xout[3] = x0[1] - x0[2]; xout[4] = x1[1] - x1[2];
    return xout;

  } else {

    // need to compute complex coefficients y0, y1.
    vector[M - 1] y0;
    vector[M - 1] y1;
    for (m in 1 : (M - 1)) {
      y0[m] = y[m] - y[2*M - m];
      y1[m] = - (x[m + 1] - x[2*M - m + 1]) * sin(2*pi() * m / N) + (y[m] + y[2*M - m]) * cos(2*pi() * m / N);
    }
    return interleave(2*M, fft_sym(2*M, x0, y0), fft_sym(2*M, x1, y1));
  }
}

vector matern_eigen(int M_nonzero, int N, real nu, real lengthscale, real delta) {
  // M_nonzero: number of non-zero frequencies
  // N: number of freqs total

  vector[M_nonzero] out;
  for (j in 1:M_nonzero) {
    out[j] = (2*nu + (lengthscale/delta)^2 * sin(pi() * j / N)^2)^(-(nu + .5) / 2);
  }
  return out; //sqrt(tgamma(nu) / tgamma(nu+.5)) * out;  //  (2*nu)^(nu/2) *
  // ((lengthscale/delta) ^ nu) * sqrt(nu) *
}

vector transform_to_matern(vector re, vector im, int N, real nu, real lengthscale, real bin_width, real shift) {
  // transforms non-centered (frequency-domain) representation with zero DC
  // component into the time-domain GP.
  // Use the 'shift' parameter
  vector[N] x;
  vector[num_elements(re)] eigens;
  eigens = matern_eigen(num_elements(re), N, nu, lengthscale, bin_width);
  x = ifft_shift_centered(N / 2 - 1, re .* eigens, im .* eigens, shift / bin_width);
  return x;
}


double
__powidf2 (double x, int m)
{
  unsigned int n = m < 0 ? -m : m;
  double y = n % 2 ? x : 1;
  while (n >>= 1)
    {
      x = x * x;
      if (n % 2)
        y = y * x;
    }
  return m < 0 ? 1/y : y;
}


"""Defines the @memoize decorator."""
from collections.abc import Hashable
from functools import partial

class memoized(object):
   """Decorator. Caches a function's return value each time it is called.
   If called later with the same arguments, the cached value is returned
   (not reevaluated).
   """
   def __init__(self, func):
      self.func = func
      self.cache = {}

   def __call__(self, *args):
      if not isinstance(args, Hashable):
         # uncacheable. a list, for instance.
         # better to not cache than blow up.
         return self.func(*args)
      if args in self.cache:
         return self.cache[args]
      else:
         value = self.func(*args)
         self.cache[args] = value
         return value

   def __repr__(self):
      """Return the function's docstring."""
      return self.func.__doc__
      
   def __get__(self, obj):
      """Support instance methods."""
      return partial(self.__call__, obj)

# Code from Stackoverflow
class Memoize:
    def __init__(self, f):
        self.f = f
        self.memo = {}

    def __call__(self, *args):
        return self.memo.setdefault(args, self.f(*args))
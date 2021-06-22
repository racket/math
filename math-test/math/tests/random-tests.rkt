#lang typed/racket
(require typed/rackunit)

(check-true
  (for/or ([n (in-range 100)])
    (> (random-natural (expt 2 64)) 4294967087)))

(for ([n (in-range 100)])
  (define x (random-natural 1 100))
  (check-true (and (<= 1 x) (< x 100))))

(for ([n (in-range 100)])
  (define x (random-bits n))
  (check-true (<= 0 x) (< x (expt 2 n))))

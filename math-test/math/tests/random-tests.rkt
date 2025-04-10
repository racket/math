#lang typed/racket
(require math/base)
(require typed/rackunit)

(define g (current-pseudo-random-generator))

(check-true
  (for/or ([n (in-range 100)])
    (> (random-natural (expt 2 64)) 4294967087)))

(for ([n (in-range 100)])
  (define x (random-integer 1 100))
  (check-true (and (<= 1 x) (< x 100))))

(for ([n (in-range 100)])
  (define x (random-bits n))
  (check-true (and (<= 0 x) (< x (expt 2 n)))))


(check-true
  (for/or ([n (in-range 100)])
    (> (random-natural (expt 2 64) g) 4294967087)))

(for ([n (in-range 100)])
  (define x (random-integer 1 100 g))
  (check-true (and (<= 1 x) (< x 100))))

(for ([n (in-range 100)])
  (define x (random-bits n g))
  (check-true (and (<= 0 x) (< x (expt 2 n)))))

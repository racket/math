#lang typed/racket/base

(require math/base
         math/array
         (only-in math/matrix matrix build-matrix)
         math/matrix/algebra
         "test-utils.rkt")

(: random-matrix (case-> (Integer Integer -> (Matrix Integer))
                         (Integer Integer Integer -> (Matrix Integer))
                         (Integer Integer Integer Integer -> (Matrix Integer))))
;; Generates a random matrix with Natural elements < k. Useful to test properties.
(define random-matrix
  (case-lambda
    [(m n)  (random-matrix m n 100)]
    [(m n k)  (array-strict (build-matrix m n (λ (i j) (random-natural k))))]
    [(m n k0 k1)  (array-strict (build-matrix m n (λ (i j) (random-integer k0 k1))))]))

(: empty-matrix (Matrix Nothing))
(define empty-matrix
  (build-simple-array #(0 0) (λ (js)
                               (error "this procedure should never be called"))))

;; ===================================================================================================
;; Determinant

(: matrix-determinant ((Matrix Exact-Rational) -> Exact-Rational))
(: matrix-determinant/row-reduction ((Matrix Exact-Rational) -> Exact-Rational))
(define-values (matrix-determinant matrix-determinant/row-reduction)
  (let ([+F (ann + (Exact-Rational * -> Exact-Rational))]
        [-F (ann - (Exact-Rational -> Exact-Rational))]
        [*F (ann * (Exact-Rational * -> Exact-Rational))]
        [/F (ann / (Exact-Rational -> Exact-Rational))]
        [=F (ann = (Exact-Rational Exact-Rational -> Boolean))])
    (: ~F (Exact-Rational Exact-Rational -> Boolean))
    (define (~F a b) (> (abs a) (abs b)))
    (values (make-matrix-determinant +F -F *F /F =F ~F)
            (make-matrix-determinant/row-reduction +F -F *F /F =F ~F))))

(check-equal? (matrix-determinant (matrix [[3]])) 3)
(check-equal? (matrix-determinant (matrix [[3]])) 3)
(check-equal? (matrix-determinant (matrix [[1 2] [3 4]])) (- (* 1 4) (* 2 3)))
(check-equal? (matrix-determinant (matrix [[1 2 3] [4  5 6] [7 8 9]])) 0)
(check-equal? (matrix-determinant (matrix [[1 2 3] [4 -5 6] [7 8 9]])) 120)
(check-equal? (matrix-determinant (matrix [[1 2 3 4]
                                           [-5 6 7 8]
                                           [9 10 -11 12]
                                           [13 14 15 16]]))
              5280)

;; TODO Determinant of the empty matrix
;; (check-equal? (matrix-determinant empty-matrix) 1)
;; (check-equal? (matrix-determinant/row-reduction empty-matrix) 1)
(for ([_  (in-range 100)])
  (define a (random-matrix 3 3 -3 4))
  (check-equal? (matrix-determinant/row-reduction a)
                (matrix-determinant a)))

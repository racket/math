#lang typed/racket/base

(require racket/list
         "matrix-types.rkt"
         "matrix-constructors.rkt"
         "matrix-basic.rkt"
         "matrix-gauss-elim.rkt"
         "utils.rkt"
         "algebra/matrix-solve.rkt"
         "../array/array-indexing.rkt"
         "../array/array-struct.rkt")

(provide
 matrix-determinant
 matrix-determinant/row-reduction  ; for testing
 matrix-invertible?
 matrix-inverse
 matrix-solve)

;; ===================================================================================================
;; Determinant

(: matrix-determinant (case-> ((Matrix Flonum) -> Flonum)
                              ((Matrix Real) -> Real)
                              ((Matrix Float-Complex) -> Float-Complex)
                              ((Matrix Number) -> Number)))
(define (matrix-determinant M)
  (define x00 (matrix-ref M 0 0))
  (define +F (add x00))
  (define -F (sub x00))
  (define *F (mul x00))
  (define /F (div x00))
  (define =F (eqv x00))
  ((make-matrix-determinant +F -F *F /F =F magnitude-compare) M))

(: matrix-determinant/row-reduction (case-> ((Matrix Flonum) -> Flonum)
                                            ((Matrix Real) -> Real)
                                            ((Matrix Float-Complex) -> Float-Complex)
                                            ((Matrix Number) -> Number)))
(define (matrix-determinant/row-reduction M)
  (define x00 (matrix-ref M 0 0))
  (define +F (add x00))
  (define -F (sub x00))
  (define *F (mul x00))
  (define /F (div x00))
  (define =F (eqv x00))
  ((make-matrix-determinant/row-reduction +F -F *F /F =F magnitude-compare) M))

;; ===================================================================================================
;; Inversion

(: matrix-invertible? ((Matrix Number) -> Boolean))
(define (matrix-invertible? M)
  (and (square-matrix? M)
       (not (zero? (matrix-determinant M)))))

(: matrix-inverse (All (A) (case-> ((Matrix Flonum)        -> (Matrix Flonum))
                                   ((Matrix Flonum) (-> A) -> (U A (Matrix Flonum)))
                                   ((Matrix Real)        -> (Matrix Real))
                                   ((Matrix Real) (-> A) -> (U A (Matrix Real)))
                                   ((Matrix Float-Complex)        -> (Matrix Float-Complex))
                                   ((Matrix Float-Complex) (-> A) -> (U A (Matrix Float-Complex)))
                                   ((Matrix Number)        -> (Matrix Number))
                                   ((Matrix Number) (-> A) -> (U A (Matrix Number))))))
(define matrix-inverse
  (case-lambda
    [(M)  (matrix-inverse M (λ () (raise-argument-error 'matrix-inverse "matrix-invertible?" M)))]
    [(M fail)
     (define m (square-matrix-size M))
     (define x00 (matrix-ref M 0 0))
     (define I (identity-matrix m (one* x00) (zero* x00)))
     (define-values (IM^-1 wps) (parameterize ([array-strictness #f])
                                  (matrix-gauss-elim (matrix-augment (list M I)) #t #t)))
     (cond [(and (not (empty? wps)) (= (first wps) m))
            (submatrix IM^-1 (::) (:: m #f))]
           [else  (fail)])]))

;; ===================================================================================================
;; Solving linear systems

(: matrix-solve
   (All (A) (case->
             ((Matrix Flonum) (Matrix Flonum)        -> (Matrix Flonum))
             ((Matrix Flonum) (Matrix Flonum) (-> A) -> (U A (Matrix Flonum)))
             ((Matrix Real) (Matrix Real)        -> (Matrix Real))
             ((Matrix Real) (Matrix Real) (-> A) -> (U A (Matrix Real)))
             ((Matrix Float-Complex) (Matrix Float-Complex)        -> (Matrix Float-Complex))
             ((Matrix Float-Complex) (Matrix Float-Complex) (-> A) -> (U A (Matrix Float-Complex)))
             ((Matrix Number) (Matrix Number)        -> (Matrix Number))
             ((Matrix Number) (Matrix Number) (-> A) -> (U A (Matrix Number))))))
(define matrix-solve
  (case-lambda
    [(M B)  (matrix-solve M B (λ () (raise-argument-error 'matrix-solve "matrix-invertible?" 0 M B)))]
    [(M B fail)
     (define m (square-matrix-size M))
     (define-values (s t) (matrix-shape B))
     (cond [(= m s)
            (define-values (IX wps) (parameterize ([array-strictness #f])
                                      (matrix-gauss-elim (matrix-augment (list M B)) #t #t)))
            (cond [(and (not (empty? wps)) (= (first wps) m))
                   (submatrix IX (::) (:: m #f))]
                  [else  (fail)])]
           [else
            (error 'matrix-solve
                   "matrices must have the same number of rows; given ~e and ~e"
                   M B)])]))

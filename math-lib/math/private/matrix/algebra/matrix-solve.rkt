#lang typed/racket/base

(require racket/fixnum
         racket/match
         "utils.rkt"
         "../matrix-types.rkt"
         "../matrix-conversion.rkt"
         "../matrix-basic.rkt"
         "../../vector/vector-mutate.rkt"
         "../../array/mutable-array.rkt")

(provide
 make-matrix-determinant
 make-matrix-determinant/row-reduction  ; for testing
 )

;; ===================================================================================================
;; Determinant

(: make-matrix-determinant
   (All (A) ((A * -> A) (A -> A)
             (A * -> A) (A -> A)
             (A A -> Boolean)
             (A A -> Boolean)
             ->
             ((Matrix A) -> A))))
(define (make-matrix-determinant +F -F *F /F =F ~F)
  (define 1F (*F))
  (define matrix-determinant/row-reduction
    (make-matrix-determinant/row-reduction +F -F *F /F =F ~F))

  (: matrix-determinant ((Matrix A) -> A))
  (define (matrix-determinant M)
    (define m (square-matrix-size M))
    (cond
      [(= m 0)  1F]
      [(= m 1)  (matrix-ref M 0 0)]
      [(= m 2)  (match-define (vector a b c d)
                  (mutable-array-data (array->mutable-array M)))
                (+F (*F a d) (-F (*F b c)))]
      [(= m 3)  (match-define (vector a b c d e f g h i)
                  (mutable-array-data (array->mutable-array M)))
                (+F (*F     a  (+F (*F e i) (-F (*F f h))))
                    (*F (-F b) (+F (*F d i) (-F (*F f g))))
                    (*F     c  (+F (*F d h) (-F (*F e g)))))]
      [else
       (matrix-determinant/row-reduction M)]))
  matrix-determinant)

(: make-matrix-determinant/row-reduction
   (All (A) ((A * -> A) (A -> A)
             (A * -> A) (A -> A)
             (A A -> Boolean)
             (A A -> Boolean)
             ->
             ((Matrix A) -> A))))
(define (make-matrix-determinant/row-reduction +F -F *F /F =F ~F)
  (define 0F (+F))
  (define 1F (*F))

  (: elim-rows! ((Vectorof (Vectorof A)) Index Index Index A Nonnegative-Fixnum -> Void))
  (define elim-rows! (make-elim-rows! +F -F *F /F =F))

  (: matrix-determinant/row-reduction ((Matrix A) -> A))
  (define (matrix-determinant/row-reduction M)
    (define m (square-matrix-size M))
    (define rows (matrix->vector* M))
    (let loop ([i : Nonnegative-Fixnum 0] [sign 1F])
      (cond
        [(fx< i m)
         (define-values (p pivot) (generic-find-partial-pivot rows m i i ~F))
         (if (=F pivot 0F)
             0F  ; no pivot means non-invertible matrix
             (let ([sign  (if (= i p) sign (begin (vector-swap! rows i p)  ; swapping negates sign
                                                  (-F sign)))])
               (elim-rows! rows m i i pivot (fx+ i 1))  ; adding scaled rows doesn't change it
               (loop (fx+ i 1) sign)))]
        [else
         (for/fold ([prod : A sign]) ([i (in-range m)])
           (*F prod (unsafe-vector2d-ref rows i i)))])))
  matrix-determinant/row-reduction)

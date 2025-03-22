#lang typed/racket/base

(require racket/performance-hint
         racket/string
         "../matrix-types.rkt"
         "../../unsafe.rkt"
         "../../array/array-struct.rkt")

(provide (all-defined-out))

(: format-matrices/error ((Listof (Array Any)) -> String))
(define (format-matrices/error as)
  (string-join (map (λ: ([a : (Array Any)]) (format "~e" a)) as)))

(: matrix-shapes (Symbol (Matrix Any) (Matrix Any) * -> (Values Index Index)))
(define (matrix-shapes name arr . brrs)
  (define-values (m n) (matrix-shape arr))
  (unless (andmap (λ: ([brr : (Matrix Any)])
                    (define-values (bm bn) (matrix-shape brr))
                    (and (= bm m) (= bn n)))
                  brrs)
    (error name
           "matrices must have the same shape; given ~a"
           (format-matrices/error (cons arr brrs))))
  (values m n))

(: matrix-multiply-shape ((Matrix Any) (Matrix Any) -> (Values Index Index Index)))
(define (matrix-multiply-shape arr brr)
  (define-values (ad0 ad1) (matrix-shape arr))
  (define-values (bd0 bd1) (matrix-shape brr))
  (unless (= ad1 bd0)
    (error 'matrix-multiply
           "1st argument column size and 2nd argument row size are not equal; given ~e and ~e"
           arr brr))
  (values ad0 ad1 bd1))

(: ensure-matrix (All (A) Symbol (Array A) -> (Array A)))
(define (ensure-matrix name a)
  (if (matrix? a) a (raise-argument-error name "matrix?" a)))

(: ensure-row-matrix (All (A) Symbol (Array A) -> (Array A)))
(define (ensure-row-matrix name a)
  (if (row-matrix? a) a (raise-argument-error name "row-matrix?" a)))

(: ensure-col-matrix (All (A) Symbol (Array A) -> (Array A)))
(define (ensure-col-matrix name a)
  (if (col-matrix? a) a (raise-argument-error name "col-matrix?" a)))


(: sort/key (All (A B) (case-> ((Listof A) (B B -> Boolean) (A -> B) -> (Listof A))
                               ((Listof A) (B B -> Boolean) (A -> B) Boolean -> (Listof A)))))
;; Sometimes necessary because TR can't do inference with keyword arguments yet
(define (sort/key lst lt? key [cache-keys? #f])
  ((inst sort A B) lst lt? #:key key #:cache-keys? cache-keys?))

(: unsafe-vector2d-ref (All (A) ((Vectorof (Vectorof A)) Integer Integer -> A)))
(define (unsafe-vector2d-ref vss i j)
  (unsafe-vector-ref (unsafe-vector-ref vss i) j))

;; Note: this accepts +nan.0
(define nonnegative?
  (λ: ([x : Real]) (not (x . < . 0))))

(define number-rational?
  (λ: ([z : Number])
    (cond [(real? z)  (rational? z)]
          [else  (and (rational? (real-part z))
                      (rational? (imag-part z)))])))

(begin-encourage-inline

  (: call/ns (All (A) ((-> (Matrix A)) -> (Matrix A))))
  (define (call/ns thnk)
    (array-default-strict
     (parameterize ([array-strictness #f])
       (thnk))))

  )  ; begin-encourage-inline

(: make-thread-local-box (All (A) (A -> (-> (Boxof A)))))
(define (make-thread-local-box contents)
  (let: ([val : (Thread-Cellof (U #f (Boxof A))) (make-thread-cell #f)])
    (λ () (or (thread-cell-ref val)
              (let: ([v : (Boxof A)  (box contents)])
                (thread-cell-set! val v)
                v)))))

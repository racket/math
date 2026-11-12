#lang typed/racket/base

(require racket/performance-hint
         racket/promise
         "../../flonum.rkt"
         "../unsafe.rkt"
         "../functions/incomplete-beta.rkt"
         "../functions/stirling-error.rkt"
         "impl/gamma-random.rkt"
         "impl/poisson-random.rkt"
         "normal-dist.rkt"
         "dist-struct.rkt"
         "utils.rkt")

(provide flnegative-binomial-pdf
         flnegative-binomial-cdf
         flnegative-binomial-inv-cdf
         flnegative-binomial-sample
         Negative-Binomial-Dist negative-binomial-dist
         negative-binomial-dist-size negative-binomial-dist-prob)

(define 2pi 6.283185307179586)

(: bad-params? (Flonum Flonum -> Boolean))
(define (bad-params? r p)
  (not (and (r . fl> . 0.0) (r . fl< . +inf.0)
            (p . fl>= . 0.0) (p . fl<= . 1.0))))

(: flnegative-binomial-log-pdf (Flonum Flonum Flonum -> Flonum))
(define (flnegative-binomial-log-pdf r p k)
  (cond [(or (bad-params? r p) (not (flinteger? k)))  +nan.0]
        [(k . fl< . 0.0)  -inf.0]
        [(fl= p 0.0)  -inf.0]
        [(fl= p 1.0)  (if (fl= k 0.0) 0.0 -inf.0)]
        [(fl= k 0.0)  (fl* r (fllog p))]
        [else
         (define n (fl+ r k))
         (define d (fl- (fl* k p) (fl* r (fl- 1.0 p))))
         (fl+ (fl+ (fl+ (fllog (fl/ r n))
                        (fl- (fl- (flstirling n) (flstirling r)) (flstirling k)))
                   (fl* 0.5 (fllog (fl/ n (fl* (fl* 2pi r) k)))))
              (fl+ (fl* r (fllog1p (fl/ d r)))
                   (fl* k (fllog1p (fl/ (fl* -1.0 d) k)))))]))

(: flnegative-binomial-pdf (Flonum Flonum Flonum Any -> Flonum))
(define (flnegative-binomial-pdf r p k log?)
  (define lp (flnegative-binomial-log-pdf r p k))
  (if log? lp (flexp lp)))

(: flnegative-binomial-cdf (Flonum Flonum Flonum Any Any -> Flonum))
(define (flnegative-binomial-cdf r p k log? 1-p?)
  (cond [(bad-params? r p)  +nan.0]
        [(k . fl< . 0.0)  (flprobability 0.0 log? 1-p?)]
        [(fl= p 0.0)  (flprobability (if (fl= k +inf.0) 1.0 0.0) log? 1-p?)]
        [(fl= p 1.0)  (flprobability 1.0 log? 1-p?)]
        [(fl= k +inf.0)  (flprobability 1.0 log? 1-p?)]
        [else
         (let ([k  (flfloor k)]
               [upper?  (if 1-p? #t #f)])
           (cond [log?  (fllog-beta-inc r (fl+ k 1.0) p upper? #t)]
                 [else  (flbeta-inc r (fl+ k 1.0) p upper? #t)]))]))

(: flnegative-binomial-inv-cdf (Flonum Flonum Flonum Any Any -> Flonum))
(define (flnegative-binomial-inv-cdf r p pr log? 1-p?)
  (cond [(or (bad-params? r p) (not (flprobability? pr log?)))  +nan.0]
        [(fl= p 1.0)  0.0]
        [(fl= p 0.0)  (if (flprobability-zero? pr log? 1-p?) 0.0 +inf.0)]
        [(flprobability-one? pr log? 1-p?)  +inf.0]
        [(flprobability-zero? pr log? 1-p?)  0.0]
        [else
         (define mean (fl/ (fl* r (fl- 1.0 p)) p))
         (define sd (fl/ (flsqrt (fl* r (fl- 1.0 p))) p))
         (define start (flmax 0.0 (flnormal-inv-cdf mean sd pr log? 1-p?)))
         (cond [1-p?  (flfind-least-integer
                       (λ: ([k : Flonum]) ((flnegative-binomial-cdf r p k log? 1-p?) . fl< . pr))
                       0.0 +inf.0 start)]
               [else  (flfind-least-integer
                       (λ: ([k : Flonum]) ((flnegative-binomial-cdf r p k log? 1-p?) . fl>= . pr))
                       0.0 +inf.0 start)])]))

(: flnegative-binomial-sample (Flonum Flonum Integer -> FlVector))
(define (flnegative-binomial-sample r p n)
  (cond [(n . < . 0)  (raise-argument-error 'flnegative-binomial-sample "Natural" 2 r p n)]
        [(bad-params? r p)  (build-flvector n (λ (_) +nan.0))]
        [(fl= p 1.0)  (build-flvector n (λ (_) 0.0))]
        [(fl= p 0.0)  (build-flvector n (λ (_) +inf.0))]
        [else
         (define ls (flgamma-sample r (fl/ (fl- 1.0 p) p) n))
         (build-flvector
          n (λ (i) (unsafe-flvector-ref (flpoisson-sample (unsafe-flvector-ref ls i) 1) 0)))]))

(define-real-dist: negative-binomial-dist Negative-Binomial-Dist
  negative-binomial-dist-struct ([size : Flonum] [prob : Flonum]))

(begin-encourage-inline
  
  (: negative-binomial-dist (case-> (-> Negative-Binomial-Dist)
                                    (Real -> Negative-Binomial-Dist)
                                    (Real Real -> Negative-Binomial-Dist)))
  (define (negative-binomial-dist [r 1.0] [p 0.5])
    (let ([r  (fl r)] [p  (fl p)])
      (define pdf (opt-lambda: ([k : Real] [log? : Any #f])
                    (flnegative-binomial-pdf r p (fl k) log?)))
      (define cdf (opt-lambda: ([k : Real] [log? : Any #f] [1-p? : Any #f])
                    (flnegative-binomial-cdf r p (fl k) log? 1-p?)))
      (define inv-cdf (opt-lambda: ([pr : Real] [log? : Any #f] [1-p? : Any #f])
                        (flnegative-binomial-inv-cdf r p (fl pr) log? 1-p?)))
      (define sample (case-lambda:
                       [()  (unsafe-flvector-ref (flnegative-binomial-sample r p 1) 0)]
                       [([m : Integer])  (flvector->list (flnegative-binomial-sample r p m))]))
      (negative-binomial-dist-struct
       pdf sample cdf inv-cdf
       0.0 +inf.0 (delay (flnegative-binomial-inv-cdf r p 0.5 #f #f))
       r p)))
  
  )

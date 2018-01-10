(require memoize)
;; polynomial will be defined as a list of terms where 
;   degree is a nonnegative integer
;   coefficient is floating point
;   terms with 0 coefficients do not need to be included
;   terms need to be sorted with first term that of the highest degree.
;   the empty list might be equal to 0.
;; need a reduce function 
;; a sort function? 

;; define the tuple list with the car as the 
(struct term (coeff deg))

(define/memo (power x n)
             (cond
               ((not (nonnegative-integer? n))
                 (error "power for positive integer exponents only"))
               ((= n 0) 1)
               (else (* x (power x (- n 1))))))                 

(define (term-eval t x0)
  (* (term-coeff t) (power x0 (term-deg t))))

(define (evaluatePoly Px x0)
  (foldr + 0 (map term-eval Px)))
  
(define (term-sort terms)
  (sort terms (lambda (x y) (> (term-deg x) (term-deg y)))))

(define (Poly+ px qx)
  (define (PolyAdd ax bx sx)
    (cond 
      ((and (empty? ax) (empty? bx)) sx)
      ((empty? bx) (PolyAdd bx ax sx)) ;;only want to bother with base case once
      ((empty? ax) (foldl cons sx bx) ;; need to double check this.
      


(define (reduce px)
  (define (inner-reduce inPol outPol)
    (if (empty? inPol) outPol
      (let ( [deg-inPol (term-deg (car inPol))]
             [deg-outPol (term-deg (car outPol))] )
        (if (= deg-inPol deg-outPol)
          (let ([hd-term (term (+ (term-coeff (car inPol))
                                  (term-coeff (car outPol)))
                               deg-inPol)])
            (inner-reduce (cdr inPol) (cons (hd-term (cdr outPol)))))
          (inner-reduce (cdr inPol) (cons (car inPol) outPol))))))
  (reverse (inner-reduce (term-sort px) '())))

(define (evalPoly Px x0)
  (define (evalP terms d x0^d sum)
    (if (empty? terms) sum 
      (let ([x^n (* x0^d (expt x0 (- (term-deg (car terms)) d)))])
        (evalP (cdr terms) 
               (term-deg (car terms)) 
               x^n 
               (+ sum (* x^n (term-coeff (car terms))))))))
  (evalP (reverse Px) 0 1 0))


;(define construct-polynomial
;  ; check terms are nonnegative integers
;  ; sort terms by degree
;  ; reduce terms
;  ; maybe: insert terms of w/ coeff 0 for missing term degrees.
;  (define (correct-terms termlist deg polyaccum)
;    (if (empty? termlist) polyaccum
;      (let* ([headterm (car termlist)]
;             [tailterms (cdr termlist)]
;             [tdeg (term-deg headterm)]
;             )
;        (case 
;          [(not (exact-nonnegative-integer? tdeg))
;           (error "not a valid polynomial")]
;          [(< tdeg deg) ;;there are two terms of the same degree.
;           (let ([polyhd (car polyaccum)]
;                 [polytail (cdr polyaccum)])
;             (construct-polynomial tailterms deg
;                  (cons (term-add polyhd headterm) polytail)))]
;          [
;

(struct polynomial (ls)
  
  (define/match poly-type?
    [(list (? pairs?) ...) 'pairpoly]
    [(list (? number) ...) 'listpoly]
    [error "invalid input"])
  )  

(struct tree (val left right))


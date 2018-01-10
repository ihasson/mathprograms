(require math)
(require infix)

(define (bisec-method f a0 b0 ε n)
  (define (bisec a b i)
    (let* ([p (/ (+ a b) 2)]
           [_ (displayln 
              (list "p =" p " , i =" i " , f(p) =" (f p)))])
      (if (or (< (abs (f p)) ε) (>= i n)) 
        p
        (if (> (* (f a) (f p)) 0)
              (bisec p b (+ i 1))
              (bisec a p (+ i 1))))))
  (bisec a0 b0 0))

;;there is something wrong with stopping based on tolerance
;;for now need to always set TOL to 0.0
(define (fixed-pt-iter fnct p0 TOL maxiter)  
  (define (fpi p i)
    (let* ( [pn (fnct p)]
            [_ (displayln (list i " : " p))] )
      (cond 
        [(>= i maxiter) 
            (begin (displayln 
              (list "Max iterations " i  " reached. Last value " p)) 
                   pn)]
        [(< (magnitude (- pn  p)) TOL) (begin (displayln pn i) pn)]
        [else (fpi pn (+ i 1)) ])))
  (fpi p0 0))

(define newtonMethod
  (λ (f fprime p0 tol n)
     (let ([g (λ (x) (- x (/ (f x) (fprime x))))])
       (fixed-pt-iter g p0 tol n))))

(define (secantMethod f a b ε n)
  (define (secantIter p0 p1 q0 q1 i)
    (let* ( [p (- p1 (* q1 (/ (- p1 p0) (- q1 q0))))]
            [_ (displayln (list i " : " p))] )
      (if (or (< (abs (- p0 p1)) ε) (> i n))
        p
        (secantIter p1 p q1 (f p) (+ i 1)))))
  (secantIter (fl a) (fl b) (f (fl a)) (f (fl b)) 2))

(define (falsePosition f a b ε n)
  (define (fpiter p0 p1 q1 i)
    (let* ( [p (- p1 (* q1 (/ (- p1 p0) (- q1 q0))))]
            [_ (displayln (list i " : " p))] )
      (if (or (< (abs (- p0 p1)) ε) (> i n))
        p
        (let ([q (f q)])
          (if (negative? (* q q1))
            (fpiter p1 p q1 q (+ i 1))
            (fpiter p0 p q0 q (+ i 1)))))))
  (fpiter (fl a) (fl b) (f (fl a)) (f (fl b)) 2))

;; 
;; differentiation 
;; taken from SICP
(define (variable? x) (symbol? x))

(define (same-variable? v1 v2)
  (and (variable? v1) (variable? v2) (eq? v1 v2)))

(define (make-sum a1 a2) (list '+ a1 a2))

(define (make-product m1 m2) (list '* m1 m2))

(define (sum? x)
  (and (pair? x) (eq? (car x) '+)))

(define (addend s) (cadr s))

(define (augend s) (caddr s))

(define (product? x)
  (and (pair? x) (eq? (car x) '*)))

(define (multiplier p) (cadr p))

(define (multiplicand p) (caddr p))

(define (deriv expr var)
  (cond ((number? expr) 0)
        ((variable? expr)
         (if (same-variable? expr var) 1 0))
        ((sum? expr)
         (make-sum (deriv (addend expr) var)
                   (deriv (augend expr) var)))
        ((product? expr)
         (make-sum
           (make-product (multiplier expr)
                         (deriv (multiplicand expr) var))
           (make-product (deriv (multiplier expr) var)
                         (multiplicand expr))))
        (else
          (error "unknown expression type -- DERIV" expr))))
;;
;; my code again!
(define (deriv-at expr var p)
  (let* ([df/dx (deriv expr var)]
         [l '(lambda (x))]
         [df (append l (cons df/dx '()))])
    (apply (eval df) (list p))))
(define (derivative expr)
  (let* ([df/dx (deriv expr 'x )]
         [l '(lambda (x))]
         [df (append l (cons df/dx '()))])
    (lambda (x) (apply (eval df) (list x)))))
   
(define g '(+ (* x (* x x)) (* -2 (* x x)) -5))  

(define (modifiedNewton f p0 ε n)
  [define fx (lambda (x) (apply (eval 
                                  (append '(lambda (x)) 
                                          f )) 
                                (list x)))]
  [define df (derivative f)]
  [define ddf (derivative (deriv f 'x))]
  [define gx (lambda (x) (- x 
                            (/ (* (df x) (fx x)) 
                               (- (* (df x) (df x)) 
                                  (* (fx x) (ddf x))))))]
  (begin (display (list (fx p0) (df p0) (ddf p0) (gx p0)))
    (fixed-pt-iter gx p0 ε n)))


;; Evaluates the polynomial P(x) and its derivative at x0.
;; returns the values P(x0) and P'(x0) not the polynomials.
(define (hornersMethod Px x0)
  (define (hm y z pxs)  
    (if (= (length pxs) 1)
      (let* ([yn (+ (* x0 y) (car pxs))]
            [_ (displayln (list "P(" x0 ") = " yn " P'(" x0 ") = " z))])
        (list yn  z))
      (let* ([yn (+ (* x0 y) (car pxs))]
             ;[_ (displayln (list yn)) ]
             )
        (hm yn (+ (* x0 z) yn) (cdr pxs)))))
  (hm (car Px) (car Px) (cdr Px)))

;;
;; passes the function x - (P(x)/P'(x)) to the fixed pt iteration
;; P(x) and P'(x) get evaluated at x by Horner's method.
;; 
(define (hornerNewton Px x0 n)
  (define g 
    (lambda (x)  
      (let ([pq (hornersMethod Px x)])
        (- x (/ (car pq) (cadr pq))))))
  (fixed-pt-iter g x0 0.0 n))

;; Divides a polynomial 
(define (synthDiv Px x0)
 ; '( an an-1           an-2 ...)
 ; '( 0  x*an 
 ; '( an a'n-1'+(x*an)    
  (define (synthD px bn Qx)
    (if (empty? px) 
      (reverse Qx)
      (let ([bn-1 (+ (car px) bn)])
        (synthD (cdr px) (* x0 bn-1) (cons bn-1 Qx)))))
  (synthD Px 0 '()))

;; For a moment I thought my Horner's Method was wrong and wrote this
;; for debugging purposes.
;; It turned that the algorithm for determining P(x) and P'(x)
;; at x are distinct from polynomial division.
(define (make-poly cs)
  (define (poly ls xn x)
    (cond 
      [(empty? ls) (error "undefined polynomial")]
      [(= (length ls) 1) (car ls)]
      [else 
        (+ (* xn (car ls)) (poly (cdr ls) (* xn x) x))]))
  (lambda (p) (poly (reverse cs) 1 p)))    

(define (roots a b c)
  (let ([part (sqrt ($ "b^2 - 4*a*c"))])
    (list ($ "(-b + part)/(2*a)") ($ "(-b - part)/(2*a)"))))

(define (MuellersMethod f p0 p1 p2 ε n)
  (define (mueller p0 p1 p2 h1 h2 δ1 δ2 d i)
    (if (> i n) 
      (begin 
        (displayln (list "method failed after " n " iterations."))
          p2)
      (let* ([b (+ δ2 (* h2 d))]
             [D (sqrt (- (* b b) (* 4 (f p2) d)))]
             [E   (if (< (magnitude (- b D)) (magnitude (+ b D)))
                    (+ b D)
                    (- b D))]
             [h (/ (* -2 (f p2)) E)]
             [p (+ p2 h)])
             
        (if (< (magnitude h) ε) 
          p
          (let* ([p0* p1]
                 [p1* p2]
                 [p2* p]
                 [h1* (- p1* p0*)]
                 [h2* (- p2* p1*)]
                 [δ1* (/ (- (f p1*) (f p0*)) h1*)]
                 [δ2* (/ (- (f p2*) (f p1*)) h2*)]
                 [d* (/ (- δ2* δ1*) (+ h2* h1*))]
                 [i* (+ i 1)])
              (mueller p0* p1* p2* h1* h2* δ1* δ2* d* i*))))))
  (let*   ([h1 (- p1 p0)]
           [h2 (- p2 p1)]
           [δ1 (/ (- (f p1) (f p0)) h1)]
           [δ2 (/ (- (f p2) (f p1)) h2)]
           [d (/ (- δ2 δ1) (+ h2 h1))]
           [i 3]) 
    (mueller p0 p1 p2 h1 h2 δ1 δ2 d i)))


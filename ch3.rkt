(require math)
(require infix)
(require plot)
(require memoize)

(plot-new-window? #t)


(define (logb_x_ b x)
  (/ (log x) (log b)))

;;Lagrange's interpolating polynomial of 1st degree
(define (P1x f x xa xb)
  (+
    (* (f xa) (/ (- x xb) (- xa xb)))
    (* (f xb) (/ (- x xa) (- xb xa)))))

;;Lagrange's interpolating polynomial of 2nd degree
(define (P2x f x xa xb xc)
  (+ 
    (* (f xa) [/ (*(- x  xb) (- x xc)) (* (- xa xb) (- xa xc))])
    (* (f xb) [/ (*(- x  xa) (- x xc)) (* (- xb xa) (- xb xc))])
    (* (f xc) [/ (*(- x  xa) (- x xb)) (* (- xc xa) (- xc xb))])))

;(define (f_2b x) (expt (- x 1) 1/3))
;(define (f_2c x) (logb_x_ 10 (- (* 3 x) 1)))
;(define (f_2d x) (- (exp (* 2 x)) x))
;(plot (function (lambda (x) (sin (* pi x))) 0 1))

;Lagrange Interp
; Rudimentry implementation of of Lagrange's interpolating polynomial. 
; Use the function Lagrange-poly to make the polynomial.
; polyfun makes the Lagrange polynomial into a function
; as an example:
; Let P be the returned value of Lagrange-poly.
; Then ((polyfun P) num) will be the lagrange poly at num.
;
(define (make-term m1 m2) (list '- m1 m2))

(define (poly-from-roots ls)
  (cons '* (map (lambda (a) (make-term 'x a)) ls)))

; Doesn't work on the Lagrange-poly directly.
; Need to figure out why.
(define (polynomial->function var poly)
  (eval (list 'lambda (cons var '())  poly)))

(define (topL v ls)
  (poly-from-roots [filter-map (λ (x) (and (not (= v x)) x)) ls]))

(define (bottomL v ls)
  ((polynomial->function 
     'x 
     (poly-from-roots [filter-map (λ (x) (and (not (= v x)) x)) ls]))
   v))

(define (Lagrange-poly xs f)
  (cons '+ (map (λ (v) 
          (list '* 
                (f v) 
                (list '/ 
                      (topL v xs) 
                      (bottomL v xs)))) 
       xs)))
; P[x] -> FN(x)
(define (polyfun Px)
  (eval (list 'lambda '(x) Px)))


;(define (f_5a x)
;  (define h
;    (dict-set* #hash()
;                8.1  16.94410
;                8.3  17.56492
;                8.6  18.50515
;                8.7  18.82091
;              ))
;  (dict-ref h x))
;
;(define p5a_P1 (Lagrange-poly '(8.3 8.6) f_5a))
;(define p5a_P2 (Lagrange-poly '(8.1 8.3 8.6) f_5a))
;(define p5a_P3 (Lagrange-poly '(8.1 8.3 8.6 8.7) f_5a))
;> ((polyfun p5a_P1) 8.4)
;17.87833
;> ((polyfun p5a_P2) 8.4)
;17.877129999999998
;> ((polyfun p5a_P3) 8.4)
;17.8771425
;

;(define (f_5b x)
;  (define h
;    (dict-set* #hash()
;               -0.75  -0.07181250
;               -0.5   -0.02475000
;               -0.25   0.33493750
;                0.0    1.10100000
;               ))
;  (dict-ref h x))
;
;; want f(-1/3)
;(define p5b_P1 (Lagrange-poly '(-.5 -.25) f_5b))
;(define p5b_P2 (Lagrange-poly '(-.5 -.25 0.0) f_5b))
;(define p5b_P3 (Lagrange-poly '(-.75 -.5 -.25 0.0) f_5b))
;(define (p5b_actual x) ($ "x^3 + 4.001*x^2 + 4.002*x + 1.101"))
;
;(define (f_5c x)
;  (define h
;    (dict-set* #hash()
;          0.1 0.62049958
;          0.2 -0.28398668
;          0.3 0.00660095
;          0.4 0.24842440
;          ))
;  (dict-ref h x))
;
;(define p5c_P1 (Lagrange-poly '(0.2 0.3) f_5c))
;(define p5c_P2 (Lagrange-poly '(0.2 0.3 0.4) f_5c))
;(define p5c_P3 (Lagrange-poly '(0.1 0.2 0.3 0.4) f_5c))
;(define (p5c_actual x) (+ 
;                         (* x (cos x)) 
;                         (* -2 (expt x 2)) 
;                         (* 3 x)
;                         -1))
;
;
;(define p13c (Lagrange-poly '(1.0 1.1 1.3 1.4) log))

; For terms in form:
; x₀,y₀ Q₀,₀=y₀
; x₁,y₁ . 
; x₂,y₂ .
; ...   .
; xₙ,yₙ .
; 

;not sure anymore that this is implemented correctly
;nope doesnt even come close.
;;recursive neville's method
; (FN x Pts) -> P[x] 
(define (recursive-neville f xs) 
  (define/memo (P xs)
    (if (= 1 (length xs))
       (f (car xs))
       (let ([a (car xs)]
             [b (car (take-right xs 1))]
             [xs-x₀ (cdr xs)]
             [xs-xₙ (drop-right xs 1)])
             (list '/ 
                   (list '- 
                         (list '* (list '- 'x a) (P xs-xₙ)) 
                         (list '* (list '- 'x b) (P xs-x₀)))
                   (- b a))
         )))
  (P xs))
;Neville's method is only supposed to give back the values at 
;particular points.
;(define Q₀ (make-hash [map (λ (xy) (map fl xy)) xy-pairs]))
;; book's version of Neville's Iterated Interpolation

;Neville-table
;  points,x -> table of values at x

;Neville's method
; points,x -> P at x
;
;(define (divided-diff-table xy-pairs)
; (make-hash [map (λ (xy) (map fl xy)) xy-pairs]))

; for each f:fn × ls(length n) → list(length n-1) 
;          f(fn , ls) ↦ (list fn(x₀,x₁),...,fn(xₙ₋₁,xₙ))
;
(define (pair-mapl fn ls)
  (cond 
    ([empty? ls] rs)
    ([= 1 (length ls)] (error "list length one"))
    ([= (length ls) 2] (list (fn (car ls) (cadr ls))))
    (else (cons (fn (car ls) (cadr ls)) 
                (pair-mapl fn (cdr ls))))))

;what is maxK? probably max number of terms or points or degree

(define (forward-diff Px ys k maxK)
  (cond 
    ([empty? ys] [error "empty list not allowed"])
    ([or (= k maxK) (= 1 (length ys))] ; the last case
       [cons (car ys) Px])
    (else 
      [forward-diff (cons (car ys) Px)
                    (pair-mapl (λ (y₀ y₁) (- y₁ y₀)) ys)
                    (+ 1 k)
                    maxK])))

(define (n-poly roots n)
  (if (= n 0) 1.0
      (cons '*
        (map (λ (x) (list '- 'x x)) (take roots n)))))

; h is ???                                       
; maxK????
; 
(define (newton-forward-divided-diff xs ys h maxK)
  (let ([k -1]
        [Pcoeffs (reverse (forward-diff '() ys 0 maxK))])
    (cons '+ 
      (map (λ (y) (begin 
                    (set! k (+ k 1))
                    (list '*
                          (/ y (* (factorial k) (expt h k))) 
                          (n-poly xs k))))
           Pcoeffs))))

(define prob3-3a (list '(0.0 0.25 0.5 0.75) '(1.0 1.64872 2.71828 4.48169)))

; cubic spline
;S(x)  
;Sⱼ(xⱼ) = aⱼ + bⱼ(x - xⱼ)¹ + cⱼ(x - xⱼ)² + dⱼ(x - xⱼ)³
;
{define (natural-cubic-spline xs f)
  [define/memo (α i)
    (if (or (<= i 0) (>= i (length xs)))
      (error "i supposed to be in {1,...,n-1}")
      (let ([aᵢ (f (list-ref xs i))]
            [aᵢ₊₁ (f (list-ref xs (+ i 1)))]
            [aᵢ₋₁ (f (list-ref xs (- i 1)))])
        (- (* (/ 3 (h i)) (- aᵢ₊₁ aᵢ)) (* (/ 3 (h (- i 1))) (- aᵢ aᵢ₋₁)))))]
  (define/memo (h i) (- (list-ref xs (+ i 1)) (list-ref xs i)))
  (define L (make-hash (list '(0 . 1.0))))
  (define Mu (make-hash (list '(0 . 0.0))))
  (define Z (make-hash (list '(0 . 0.0))))
  (define C (make-hash (list '(0 . 0.0))))
  (define A (make-hash))
  (define B (make-hash))
  (define D (make-hash))
  [define (solve-tri-diag-syst lᵢ₋₁ μᵢ₋₁ zᵢ₋₁ i n)
    (if (= n i)
        (begin
          ;(displayln (list "solve tri diag" i n))
          (hash-set! L n 1.0)
          (hash-set! Z n 0.0)
          (hash-set! C n 0.0)
          ;(displayln (list "C" C))
          )
        
        (let* (;[_ (displayln (list "solve tri diag" i))]
               [hᵢ (h i)]
               [hᵢ₋₁ (h (- i 1))]
               [αᵢ (α i)]
               [lᵢ (- (* 2 (- (list-ref xs (+ i 1)) (list-ref xs (- i 1))))
                      (* hᵢ₋₁ μᵢ₋₁))]
               [μᵢ (/ hᵢ lᵢ)]
               [zᵢ (/ (- αᵢ (* hᵢ₋₁ zᵢ₋₁)) lᵢ)])
          (begin
            (hash-set! L i lᵢ)
            (hash-set! Mu i μᵢ)
            (hash-set! Z i zᵢ)
            (solve-tri-diag-syst lᵢ μᵢ zᵢ (+ i 1) n))))]
  [define N (- (length xs) 1)]
  [define (set-jth j)
    (let* ([_ (displayln j)]
           [a (hash-ref A j)]
           [_ (displayln "c")]
           [_ (displayln (list Z Mu C))]
           [c (- (hash-ref Z j) (* (hash-ref Mu j) (hash-ref C (+ 1 j))))]
           [_ (displayln "b")]
           [b (- (/ (- (hash-ref A (+ 1 j)) (hash-ref A j)) (h j))
                 (/ (* (h j) (+ (hash-ref C (+ j 1)) (* 2 c))) 3.0))]
           [_ (displayln "d")]
           [d (/ (- (hash-ref C (+ 1 j)) c) (* 3.0 (h j)))])
      (begin
        (hash-set! C j c)
        (hash-set! B j b)
        (hash-set! D j d)
        (list (list-ref xs j) a b c d)))]
  (begin (map (λ (i) (hash-set! A i (f (list-ref xs i))))
              (build-list (length xs) values))
         (solve-tri-diag-syst 1.0 0.0 0.0 1 N)
         ;(displayln (list (build-list (length xs) values)))
         (map (λ (j) (set-jth j))
              (reverse (drop-right (build-list (length xs) values) 1))))}
(define (cubic-spline-function list-of-Sjs)
  (define xs (sort (map (λ (ls) (car ls)) list-of-Sjs) <))
  (define Sj-Table (make-hash (map (λ (sj) (list (car sj) sj)) list-of-Sjs)))
  (define (find-Sj ls x prev)
    (cond 
      [(empty? ls) prev]
      [(< x (car ls)) prev]
      [else (find-Sj (cdr ls) x (car ls))]))
  (define (Sj->function x_j a b c d)
    (lambda (x) (+ a 
                   (* b (- x x_j) ) 
                   (* c (expt (- x x_j) 2))
                   (* d (expt (- x x_j) 3)))))
  (lambda (x) 
    (letrec ([Sj (if (< x (car xs))
                (hash-ref Sj-Table (car xs))
                (hash-ref Sj-Table (find-Sj (cdr xs) x (car xs))))]
             [_ (displayln Sj)])
      
            ((apply Sj->function (car Sj)) x))))

; n+1 equally spaced points in the interval
(define (interval a b n)
  (let* ([Δ (- b a)]
         [h (/ Δ n)]
         )
    (build-list (+ n 1) (lambda (x) (+ a (* h x))))))


(define (three-point-endpoint y₀ y₁ y₂ h)
  (* (/ 1 (* 2 h)) [+ (* -3 y₀) (* 4 y₁) (* -1 y₂)]))

(define (three-point-midpoint y₀ y₁ y₂ h)
  (* (/ 1 (* 2 h)) [- y₂ y₀])) 

; note that absolute error is alread the the best distance function for 
; the job.
(define (closest-k-points x xs k distfun)
  (take 
    (sort xs (lambda (a b) (< (distfun x a) (distfun x b))))
    k))

;; Doesn't give the same answer as the book
;; not sure why and need to debug.
(define (Richardson-extrapolation f x h i)
  (if (= i 1) 
    (/ (- (f (+ x h)) (f x)) h )
    (/ 
      [- 
        (* (expt 4 (- i 1)) (Richardson-extrapolation f x (/ h 2) (- i 1)))
        (Richardson-extrapolation f x h (- i 1))]
      (- [expt 4 (- i 1)] 1))
    )
  )
; x = x₀ ;not working properly
(define (Rich1 f x h i)
  (if (= i 1)
    [/ (- (f (+ x h)) (f x)) h] ; the problem was with computing the deriv.
    [+ (Rich1 f x (/ h 2) (- i 1))
       (/ (- 
            (Rich1 f x (/ h 2) (- i 1))
            (Rich1 f x h (- i 1)))
          (- (expt 4 (- i 1)) 1))]))
;;; Gives an approximation of the derivative at x ;;no this doesn't 
;; Richardson's Extrapolation                     ;; work either.
; gives same answers as the book.
(define (Rich2 f x h i)
  (if (= i 1) 
    [/ (- (f (+ x h)) (f x)) h]
    [+ (Rich2 f x (/ h 2) (- i 1)) 
       (- (Rich2 f x (/ h 2) (- i 1))
          (Rich2 f x h (- i 1)))]))

(define (regular-deriv f x h)
  (/ (- (f (+ x h)) (f x)) h ))
(define (centered-diff f x h)
  (/ (- (f (+ x h)) (f (- x h))) 2 h))
(define (richardson-extrap f x h i approximation-function)
  (define (rich f x h i)
    (if (= i 1) 
      [approximation-function f x h]
      [+ (rich f x (/ h 2) (- i 1)) 
         (/ (- (rich f x (/ h 2) (- i 1))
               (rich f x h (- i 1)))
            (- (expt 4 (- i 1)) 1))]))
  (rich f x h i))

;; indexed starting at 1
(define (unzip ls)
  (define (unzip* l acc-odd acc-even)
    (cond
      [(empty? l) (list (reverse acc-odd) (reverse acc-even))]
      [(= (length l) 1) (list (reverse (cons (car l) acc-odd))
                              (reverse acc-even))]
      [else
        (unzip* (cddr l) (cons (car l) acc-odd) (cons (cadr l) acc-even))]))
  (unzip* ls '() '()))

;; composite Simpson's rule
;
(define (old-composite-simpson f a b n)
  (cond 
    [(odd? n) (error "n needs to be even not odd")]
    [(<= n 0) (error "n needs to be a positive integer.")]
    [else 
      (let* 
        ([h (/ (- b a) n)]
         [X (unzip (drop-right (drop (interval a b n) 1) 1))] 
         [X-even (cadr X)]
         [X-odd  (car X)]) 
        (* (/ h 3)
           (+
             (f a)             
             (* 2 (apply + (map f X-even)))
             (* 4 (apply + (map f X-odd)))
             (f b))))]))

(define (summation f x h i n)
  (cond
    [(> i n) 0] ;not sure wether or not to keep this line 
                ;or maybe make that condition an error
    [(= i n) (f x)]
    [else 
      (+ (f x) (summation f (+ x h) h (+ i 1) n))]))

(define (composite-simpson f a b n)
  (let ([h (/ (- b a) n)])
    (*
      (/ h 3)
      (+
        (f a)
        (* 2 (summation f (+ a h h) (+ h h) 1 (- (/ n 2) 1)))
        (* 4 (summation f (+ a h) (+ h h) 1 (/ n 2)))
        (f b)))))

(define (composite-trapezoidal f a b n)
  (let ([h (/ (- b a) n)])
    (* 
      (/ h 2)
      (+
        (f a)
        (* 2 (summation f (+ a h) h 1 (- n 1)))
        (f b)))))

(define/memo (Romberg k j f a b)
  (cond
    [(< k j) (error "k<j is invalid")]
    [(and (= j 1) (>= k j)) ; could reduce to just (= j 1)
     (composite-trapezoidal f a b (expt 2 (- k 1)))]
    [else 
      (+ 
        (Romberg k (- j 1) f a b)
        (/ (- 
             (Romberg k (- j 1) f a b)
             (Romberg (- k 1) (- j 1) f a b))
           (- (expt 4 (- j 1)) 1)))]))

(define (simpsons-rule f a b)
  (let ([h (/ (- b a) 2)])
    (* (/ h 3) (+ (f a) (f b) (* 4 (f (+ a h)))))))

;need the following:
;a function f, an interval [a,b], previous error, tolerance
;
;k to keep track of how many time the interval has been divided and a max-K
;size of interval h is derived.
(define (adaptive-quadrature f a-initial b-initial 
                             previous-error 
                             initial-tolerance 
                             max-divisions
                             accumulatedsum)
  ;rξ₀ = (absolute-error S(a,b) S(b-a/2,b))
  ;lξ₀ = (absolute-error S(a,b) S(a,b-a/2))
  ;ξi = (absolute-error (- f@b f@a) ξi-1)
  (define (subdivide? k S1 S2 current-tolerance)
      (and 
        (< k max-divisions)
        (< current-tolerance (absolute-error S1 S2))))
  (define (quadrature a b a+h S-prev k)
    (let ([Sl (simpsons-rule f a a+h)]
          [Sr (simpsons-rule f a+h b)]
          [current-tol (/ initial-tolerance (expt 2 k))]);double check this
      (+ 
        (if (subdivide? k Sl S-prev current-tol) 
          (quadrature a a+h (/ (+ a a+h) 2) Sl (+ k 1))
          Sl)
        (if (subdivide? k Sr S-prev current-tol)
          (quadrature a+h b (/ (+ a+h b) 2) Sr (+ k 1))
          Sr))))
  (let ([S (simpsons-rule f a-intial b-initial)])
    (quadrature a-initial b-initial (/ (+ a b) 2) S 0)))
  
;;this one doesn't work but other romberg integration works flawlessly.
(define (romberg-integration f a b ε max-n)
  (define/memo (Rkj k j)
    (if (= j 1)
      (composite-trapezoidal f a b (sub1 (expt 2 (sub1 k))))
      (+ (Rkj k (sub1 j))
         (/ (- (Rkj k (sub1 j)) (Rkj (sub1 k) (sub1 j)))
            (sub1 (expt 4 (sub1 j)))))))
  (define (iter i)
    (if (or 
          (= i max-n)
          (< (absolute-error (Rkj i i) (Rkj (sub1 i) (sub1 i)))
             ε))
      (Rkj i i)
      (iter (add1 i))))
  (iter 1))

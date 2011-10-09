open Printf

(* Ocaml program to calculate Pi

Based on python code by Nick Craig-Wood <nick@craig-wood.com>

See: http://www.craig-wood.com/nick/articles/pi-chudnovsky/ for more
info

Compile with:
  ocamlbuild -use-ocamlfind -pkgs zarith,unix pi.native

 *)


(* Compute int_of_float(pi * 10^digits) *)
let pi_chudnovsky_bs =
  let c = Z.of_int 640320 in
  let c3_over_24 = Z.(c * c * c / of_int 24) in
  let six = Z.of_int 6 and five = Z.of_int 5 and two = Z.of_int 2 in
  let ten = Z.of_int 10 in
  let c1 = (Z.of_int 13591409) in
  let c2 = (Z.of_int 545140134) in
  let tab_gen a pab = 
    let tab = Z.(pab * (c1 + c2 * a)) in
    if Z.(a land one = one) then Z.(~- tab) else tab
  in
  let rec bs a b =
    (* 
       Compute terms for binary splitting the Chudnovsky infinite series

       t(a) = +/- (13591409 + 545140134*a)
       p(a) = (6*a-5)*(2*a-1)*(6*a-1)
       b(a) = 1
       q(a) = a * a * a * c3_over_24

       returns P(a,b) Q(a,b) and T(a,b)
     *)
    if b - a = 1 then
      let a = Z.of_int a in
      if a = Z.zero then 
        Z.(one,one,c1)
      else
        let pab = Z.((six*a-five)*(two*a-one)*(six*a-one)) in
        pab, Z.(a*a*a*c3_over_24), tab_gen a pab
    else
      let m = (a+b)/2 in
      let pam,qam,tam = bs a m in
      let pmb,qmb,tmb = bs m b in
      Z.(pam*pmb, qam*qmb, qmb*tam + pam*tmb)
  in
  let digits_per_term = log10(Z.to_float c3_over_24 /. 6. /. 2. /. 6.) in
  fun digits -> 
    let n = 1 + int_of_float (float digits /. digits_per_term) in
    let p,q,t = bs 0 n in
    let one_squared = Z.pow ten (2 * digits) in
    let sqrt_c = Z.(sqrt (of_int 10005 * one_squared)) in
    Z.(q * of_int 426880 * sqrt_c / t)

let check_digits = [100,70679; 1_000,1989; 10_000, 75678; 100_000, 24646; 1_000_000, 58151; 10_000_000, 55897]

let main () = 
  Z.print (pi_chudnovsky_bs 100);
  print_newline();
  let upper_lim = 
    if Array.length Sys.argv > 1 then int_of_string Sys.argv.(1) else 9 
  in
  for log10_digits = 1 to upper_lim do
    let digits = Z.(to_int (pow (of_int 10) log10_digits)) in
    let start = Unix.gettimeofday() in
    let pi = pi_chudnovsky_bs digits in
    printf "chudnovsky_bs_z: digits=%d time=%.2f\n%!" digits (Unix.gettimeofday() -. start);
    try
      let check = List.assoc digits check_digits in
      let last_five = Z.(to_int (rem pi (of_int 100_000))) in
      if check = last_five then 
        printf "Last 5 digits %05d OK\n%!" last_five
      else
        printf "Last 5 digits %05d wrong should be %05d\n%!" last_five check
    with Not_found -> ()
  done

let () = main ()

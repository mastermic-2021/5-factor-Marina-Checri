encode(m)=fromdigits(Vec(Vecsmall(m)),128);
decode(m)=Strchr(digits(m, 128));
[n,c] = readvec("input.txt");


/* ****************************************************************************** */
/* ******************************* Explications : ******************************* */
/* ****************************************************************************** */

/* Basile cherche à déchiffrer un message c=[m^2, (m/n), m%2] transmis par Adèle.
 * Il possède la clef publique n=p*q et doit donc la factoriser pour obtenir
   explicitement p et q.

 * Pour ce faire, on peut utiliser les méthodes de Pollard vues en cours :
   * la méthode rho de Pollard (mais ici, elle est bien moins performante que la seconde :)
   * la méthode p-1 de Pollard.
  qui permettent toutes deux de retrouver la factorisation d'un nombre.
 * Elles sont codées dans la section suivante.

  * Puis, il lui suffit de déchiffrer le message c.
  * Rappel :
    * c = m^2 mod n
    * m_p = c^((p+1)/4) mod p
    * m_q = c^((q+1)/4) mod q

  * Et il suffit d'appliquer le lemme chinois pour retrouver n !
  * C'est ce que fait la fonction dechiffre(c,p,q).
*/

/* ****************************************************************************** */
/* ******************************** Fonctions : ********************************* */
/* ****************************************************************************** */

/* Méthode de chiffrement de Rabin. */
chiffre(m) = [m^2%n,kronecker(m,n),m%2];

/* Méthode de déchiffrement de Rabin. */
dechiffre(c,p,q)={
  \\n=p*q
  \\c=[m^2, (m/n), m%2]
  my(m_p, m_q, u, v);
  m_p = (Mod(c[1],n))^((p+1)/4);
  m_q = (Mod(c[1],n))^((q+1)/4);
  [u,v]= bezout(p,q);

  return (lift(n - (u * p * m_q - v * q * m_p)));
}


/* Méthode rho de pollard telle qu'illustrée en cours. */

Rho_de_Pollard(n)={
  my(x,y,pgcd);
  x = 1;
  y = 1;
  pgcd = 1;
  while(pgcd==1, x=(x^2+1)%n; y=(y^2+1)%n; y=(y^2+1)%n; pgcd = gcd(x-y,n););
  pgcd;
}


/* Méthode p-1 de Pollard vue en cours, avec une borne B friable et a fixé à 3 par défaut. */

pollard_p_moins_1(n,B=10000,a=3) = {
  my(pgcd);
  a = Mod(a,n);
  forprime(p=2,B,
    a=a^(p^logint(n,p));
    pgcd = gcd(lift(a)-1,n);
    if(pgcd>1,break())
    );
  pgcd;
};


/* Méthode p-1 de Pollard améliorée :
 * on ne s'arrête pas tant qu'on a pas réussi à factoriser,
 * (donc tant que le pgcd vaut 1),
 * et on augmente alors la borne.
 * Par défaut, a=3 (comme dans le fonction précédente).
*/

pollard_p_moins_1_ameliore(n,a=3) = {
  my(i,pgcd);
  a = Mod(a,n);
  i = 1;
  pgcd = 1;
  while(pgcd==1,
    a=a^i;
    pgcd = gcd(lift(a)-1,n); 
    i=i+1;
    );
  pgcd;
};



\\ On factorise n, et on déchiffre !

p=pollard_p_moins_1_ameliore(n);

m=dechiffre(c, p, n/p);
m=decode(m);

print(m);



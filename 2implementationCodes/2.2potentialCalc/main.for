
      IMPLICIT REAL(A-H, O-Z), INTEGER(I-N) 
      
      integer     i, j, k, tar
      real        atom(30000,3),charge(30000)
      real        ug(30000)
      real        cstr(30000,3000),sstr(30000,3000)
      real        potential,sht,longt,V0,css
      integer     nlocal
	   real        dx,dy,dz,rij,kcount
      real        g_ewald,prefactor,volume,back
      character   NullChar
      real        pi

      pi = 4.*atan(1.)
      tar = 3 ! the id of the atom of which the potential is calcd. 

	   open(1,file="param.txt") 
	   read(1,*) nlocal
      read(1,*) all
	   read(1,*) kcount
	   read(1,*) g_ewald
	   read(1,*) volume
	   close(1)
	  
	   open(1,file="atoms.txt")  
	   do i = 1,all
	  	   read(1,*) NullChar,atom(i, :),charge(i)
	   end do
      close(1)

	   open(1,file="ug.txt")
	   do i = 1,kcount
	  	   read(1,*) NullChar,ug(i)
	   end do
	   close(1)
	  
	   open(1,file="cossin.txt")
	   do i = 1,kcount
         do j=1,nlocal
            read(1,*) NullChar,cstr(i,j),sstr(i,j)
         end do
	   end do
	   close(1)
         
	   potential = 0.0
	   sht       = 0.0
	   longt     = 0.0
	   V0		    = 1.0

      ! short term 
	   do i = 1,all
         dx = atom(i,1) - atom(tar,1)
	  		dy = atom(i,2) - atom(tar,2)
	  		dz = atom(i,3) - atom(tar,3)
	  		rij = (dx**2+dy**2+dz**2)**0.5
         if (rij < 12.0.and.rij > 0.0) then
             sht = sht + charge(i)/rij*erfc(g_ewald*rij)
         end if
	   end do

	   ! long term
	   do i = 1,kcount
	  	   do j = 1,nlocal
            css = cstr(i,tar)*cstr(i,j)+sstr(i,tar)*sstr(i,j)
	  	      longt = longt + ug(i)*charge(j)*css
		   end do
	   end do

	   ! self term
	   self = -2.*g_ewald*charge(tar)/sqrt(pi)

	   ! background term
	   prefactor = pi/(g_ewald**2*volume)
	   do i = 1, nlocal
	  	   back = back + charge(i)
	   end do
	   back = back*prefactor

      print *, sht, longt, self, back, V0
	  
	   potential = potential + sht + longt + self + back - V0
	  
	   print *, potential
     
      open(1,file='result.txt')
      write(1,*) "potential = ", sht + longt + self + back
      write(1,*) "short-range term = ", sht
      write(1,*) "long-range term = ", longt
      write(1,*) "self-term = ", self
      write(1,*) "back-ground screening = ", back
      write(1,*) "target potential = ", V0
      write(1,*) "error = ", potential/V0
      close(1)

      end


#  MINDZ  v1.0.0  

​

The purpose of this markdown is to document the procedure for creating new classes in the model.  

​

##  Description  

​

To add a new class, you must first inherit it from the existing `Zoo` class which defines *generic* advected particles (only considers Lagrangian transport, without any biological reaction terms), then add the **attributes** and **biological  procedures** specific to the *new* class.  

​

##  Class  hierarchy  

​

![](https://i.imgur.com/actpgxT.png)  

​

In the example above, the `Krill` subclass inherits all of `Zoo` class attributes and procedures.  

​

##  Files  to  modify  

​

1.  `MINDZ/model/partSubClass_mod`  :  This  module  contains  the  variables  and  parameters  **definitions**,  **attributes**  and  **procedures**  (functions  and  subroutines)  associated  with  the  subclasses  derived  from  the  `Zoo`  class.  

​

2.  `MINDZ/model/partFun_mod`  :  This  module  implements  the  procedures  for  the  **generation**  of  the  *list  of  particles*,  their  **temporal  evolution**  according  to  *biological  functions*,  and  the  **creation**  of  the  *output  file*.  

​


3. `lagrangian_motion.f90` : This module implements the particle transport protocols with _parallelisation_. 
4.  `MINDZ/run/run.list`:  This  file  contains  the  values  of  the  simulation  parameters  that  are  red  at  runtime  in  a  **namelist**  format;  hence  the  model  *does  not  have*  to  be  compiled  again  for  changes  in  these  values  to  take  effect.  


​

##  Adding  a  new  class  

​

###  `partSubClass_mod.f90`  

​

1.  Open  the  `partSubClass_mod.f90`  source  file  

​

2.  Head  down  to  the  dedicated  section  for  new  classes  

​

3.  Declare  your  new  class  as  publicly  accessible:  `public  ::  newClass`  

​

4.  Explicitely  declare  that  your  new  class  is  **inherited**  from  the  `Zoo`  class:  `type,  extends(Zoo)  ::  newClass`  

​

5.  Declare  the  new  class  **variables**  below:  `<variable  type>  ::  variable  name  ![comments]`  

​

6.  Declare  the  **procedures**  specific  to  the  new  class  (*not*  already  part  of  the  `Zoo`  class):  

    -  Under  the  keyword  `contains`,  declare  your  class  bound  procedures:  `procedure,  <access  mode>  ::  procedure_name  ![comments]`  

​

    -  Close  the  class  declaration  section:  `end  type  newClass`  

​

7.  Under the  keyword  `contains`,   implement  your  procedures  (*getters*,  *setters*,  *constructor*  and  other). 

​

###  `run.list`  

​

Add a dedicated **namelist** section `&newClass_nml` with the parameters for your new class:  

```{FORTRAN}

&newClass_nml

 var_logical_1 = .false.,

 var_int_1     = 1,

 var_real_1    = .true.,

 var_char_1    = '../PATH/output.dat' /

```

​

###  `partFun_mod.f90`  

1.  Open  the  `partfun_mod`  source  file  

​

2.  Implement  your  new  class  generator  subroutine  `GeneratenewClass`:  

​

-  Follow  the  structure  of  `GenerateKrill`  subroutine  

​

-  Change  `init_krill`  with  your  new  class  constructor  

​

-  Change  the  name  of  the  namelist  section  to  your  own  `&newClass_nml`  

3. Implement  your  new  class  subroutine  `evolvenewClass` :
-  Copy the `evolveKrill` subroutine and substitute `Krill` instances with your newClass. (`evolveCalanus` is a good example )
​
###  `lagrangian_motion.f90`  
1.  Open  the  `lagrangian_motion`  source  file  

​

2.  Head down to the `do while` loop block in the `trajectory` subroutine, add your class type inside the `select type` block, this assigns the pointer array for OMP parallelization. 
 
		 	 select type(curr => mylist%currentValue())

			class  is (Zoo)

				ptrarray(i)%ptr => curr

			class  is (Krill)

				ptrarray(i)%ptr => curr

			class  is (Calanus)

				ptrarray(i)%ptr => curr
			
			class is (NewClass) 
		
				ptrarray(i)%ptr => curr

			end  select 
## Implementing your new class' procedures

The inherited class can implement internal procedures that are not present in the mother class. 
1. Open the `partSubClass_mod.f90` source file. 

2. Head down to your new class declaration `type, extends(Zoo) :: newClass`. 

3. The procedures are declared **inside** your new class, but are implemented **outside**: 

	- Declare your procedures under the `contains` statement located inside your new class.

	-  If you wish to assign a name for your procedure that's already been used by another class , you can do so by **overriding** your procedure, this will allow your procedure  to be referenced by a shared name while having  your specific new class implementation. 
 

 Example: 

     type, extends(Zoo) :: newClass
     ! Variable declaration
     contains
     procedure, public :: Procedure_name  => Procedure_name_newClass
 

- If `Procedure_name ` is called upon a `newClass` object,    `Procedure_name_newClass` will be called instead. 
4. You need to implement your procedures outside of your class declaration, just make sure you reference it by its **overidden** name: 
-	Under the `contains` statement (outside of class declaration), implement your procedure: 
 
		 	 

			      function/subroutine  Procedure_name_newClass(this)
				    class(newClass) :: this
				    ! Code here ! 
				   end function/subroutine   Procedure_name_newClass


   




				 


#  END


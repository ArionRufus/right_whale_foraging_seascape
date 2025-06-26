module link_mod

  private
  public :: link

  ! The link type is used to create a list in 
  ! which can be stored any type of object

  type link
     private

     class(*),    pointer :: value => null() ! value stored in link
     class(link), pointer :: next  => null() ! next link in list

   contains
     procedure :: getValue    ! return value of link
     procedure :: printLinks  ! prints linked list starting with this link
     procedure :: nextLink    ! return next pointer
     procedure :: setNextLink ! set next pointer

  end type link

  ! The general interface is used to call a
  ! function or subroutine when creating an 
  ! object of type x.
  ! In the current situation the object
  ! would be an object of type link

  interface link
     procedure constructor
  end interface link

contains

  function constructor(value, next)
    class(link), pointer :: constructor
    class(*)             :: value
    class(link), pointer :: next
    
    allocate(constructor)
    constructor%next => next
    
    allocate(constructor%value, source=value)
  end function constructor

  function getValue(this)
    class(link)          :: this
    class(*)   , pointer :: getValue

    getValue => this%value
  end function getValue

  function nextLink(this)
    class(link)          :: this
    class(link), pointer :: nextLink

    nextLink => this%next
  end function nextLink

  subroutine setNextLink(this, next)
    class(link)          :: this
    class(link), pointer :: next

    this%next => next
  end subroutine setNextLink

  subroutine printLink(this)
    class(link)          :: this

    select type(v => this%value)
    type is (integer)
       print *, v
    type is (real)
       print *, v
    class default
       stop 'printLink : unexpected type for link'
    end select
  end subroutine printLink

  subroutine printLinks(this)
    class(link)          :: this
    class(link), pointer :: curr

    call printLink(this)
    curr => this%next

    do while(associated(curr))
       call printLink(curr)
       curr => curr%next
    end do
  end subroutine printlinks

end module link_mod
		

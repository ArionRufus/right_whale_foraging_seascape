module abstList_mod

  use link_mod

  private

  public :: list

  !--------------------------------------------------------------------72
  type, abstract :: list
     private

     class(link), pointer :: firstLink => null() ! first link in list
     class(link), pointer :: lastLink  => null() ! last link in list
     class(link), pointer :: currLink  => null() ! list iterator

   contains
     procedure, non_overridable :: addValue     ! add class(*) to list
     procedure, non_overridable :: firstValue   ! returns value of first link in line
     procedure, non_overridable :: reset        ! reset list iterator
     procedure, non_overridable :: next         ! increment list iterator
     procedure, non_overridable :: currentValue ! get value from iterLink
     procedure, non_overridable :: moreValues   ! more values for iterator?

  end type list

  !--------------------------------------------------------------------72
  abstract interface
     subroutine printValues(this)
       import list
       class(list) :: this
     end subroutine printValues
  end interface


contains

  !--------------------------------------------------------------------72
  subroutine addvalue(this, value)
    class(list)          :: this
    class(*)             :: value
    class(link), pointer :: newlink

    if (.not. associated(this%firstLink)) then
       this%firstLink => link(value, this%firstLink)
       this%lastlink  => this%firstLink
    else
       newLink        => link(value, this%lastLink%nextLink())
       call this%lastLink%setNextLink(newLink)
       this%lastLink  => newLink
    endif
  end subroutine addvalue

  !--------------------------------------------------------------------72
  subroutine next(this)
    class(list)          :: this

    this%currLink => this%currLink%nextLink()
  end subroutine next

  !--------------------------------------------------------------------72
  subroutine reset(this)
    class(list)          :: this

    this%currLink => this%firstLink
  end subroutine reset

  !--------------------------------------------------------------------72
  function firstValue(this)
    class(list)          :: this
    class(*)   , pointer :: firstValue

    firstValue => this%firstLink%getValue()
  end function firstValue

  !--------------------------------------------------------------------72
  function currentValue(this)
    class(list)          :: this
    class(*)   , pointer :: currentValue

    currentValue => this%currLink%getValue()
  end function currentValue

  !--------------------------------------------------------------------72
  function moreValues(this)
    class(list)          :: this
    logical moreValues

    moreValues = associated(this%currLink)
  end function moreValues

end module abstList_mod





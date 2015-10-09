program print_itgl
character*100 w1,w2,w3,w4,w5,w6,w7,w8,w9
logical FileExist
integer nbase,nelectron
double precision,allocatable,dimension(:,:,:,:)::g4
double precision,allocatable,dimension(:)::oneh
double precision,allocatable,dimension(:)::ovrlp
double precision,allocatable,dimension(:,:)::oneh_full
double precision,allocatable,dimension(:,:)::ovrlp_full

!initialize aces
        call aces_init_rte
        call aces_ja_init 
!read n_base and print
        call getrec(20, "JOBARC", "NBASTOT", 1, nbase)
         write (*,*) nbase
!read n_base, n_alpha, n_beta
        FileExist=.false.
        inquire(file='out',exist=FileExist)
        if (FileExist) then
                open(unit=10,file='out',access='SEQUENTIAL')
                rewind 10
                do 
                        read (10,*) w1
                        !get n_base
!                        if (w1=="@VSCF:") then
!                                backspace 10
!                                read (10,*) w1,w2,w3,w4,w5,w6,w7,w8,w9
!                                if(w2=="There" &
!                                    .and. w3=="are" &
!                                    .and. w5=="functions" &
!                                    .and. w6=="in" &
!                                    .and. w7=="the"&
!                                    .and. w8=="AO"&
!                                    .and. w9=="basis."&
!                                ) then
!                                        backspace 10
!                                        read (10,*) w1,w1,w1,nbase
!                                        write (*,*) nbase
!                                end if
!                        end if
                        !get alpha and beta population
                        if (w1=="Alpha" .or. w1=="Beta") then
                                backspace 10
                                read (10,*) w1,w2,w3,w4
                        endif
                        if (w2=="population" &
                            .and. w3=="by" &
                            .and. w4=="irrep:"&
                           ) then
                                backspace 10
                                read (10,*) w1,w2,w3,w4,nelectron
                                write (*,*) nelectron
                                if (w1=="Beta") then
                                        exit
                                end if
                        end if
                end do
                close(10)
        end if
!read 1e integrals
        NDIM = nbase
        allocate (oneh((NDIM*(NDIM+1))/2))
        allocate (ovrlp((NDIM*(NDIM+1))/2))
        allocate (oneh_full(NDIM,NDIM))
        allocate (ovrlp_full(NDIM,NDIM))
        call Get1EInt(nbase,oneh,ovrlp,oneh_full,ovrlp_full)

!read 2e integrals
        allocate (g4(nbase,nbase,nbase,nbase))
        call zero(g4,nbase*4)
        call Get2EInt(nbase,g4)
        call aces_fin
end program print_itgl

program print_itgl
character*100 w1,w2,w3,w4,w5,w6,w7,w8,w9
logical FileExist
integer nbase,nelectron


!initialize aces
        call aces_init_rte
        call aces_ja_init 
!read n_base and print
!        call getrec(20, "JOBARC", "NBASTOT", 1, nbase)
!        write (*,*) nbase
!read n_base, n_alpha, n_beta
        FileExist=.false.
        inquire(file='out',exist=FileExist)
        if (FileExist) then
                open(unit=10,file='out',access='SEQUENTIAL')
                rewind 10
                do 
                        read (10,*) w1
                        !get n_base
                        if (w1=="@VSCF:") then
                                backspace 10
                                read (10,*) w1,w2,w3,w4,w5,w6,w7,w8,w9
                                if(w2=="There" &
                                    .and. w3=="are" &
                                    .and. w5=="functions" &
                                    .and. w6=="in" &
                                    .and. w7=="the"&
                                    .and. w8=="AO"&
                                    .and. w9=="basis."&
                                ) then
                                        backspace 10
                                        read (10,*) w1,w1,w1,nbase
                                        print *,nbase
                                end if
                        end if
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
                                print *,nelectron
                                if (w1=="Beta") then
                                        exit
                                end if
                        end if
                end do
                close(10)
        end if
!read 1e integrals
        call Get1EInt(nbase)
!read 2e integrals
        call Get2EInt(nbase)
        call aces_fin
end program print_itgl

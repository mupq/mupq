ELFNAME = $*

obj/_ELFNAME_%.o: $(CONFIG)
	@echo "  GEN     $@"
	$(Q)echo "const char _elf_name[] = \"$(ELFNAME)\";" | \
		$(CC) -x c -c -o $@ $(filter-out -g3,$(CFLAGS)) -

elf/%.elf: obj/_ELFNAME_%.elf.o $(LINKDEPS) $(CONFIG)
	@echo "  LD      $@"
	$(Q)[ -d $(@D) ] || mkdir -p $(@D)
	$(Q)$(LD) $(LDFLAGS) -o $@ $(filter %.o,$^) -Wl,--start-group $(LDLIBS) -Wl,--end-group

obj/%.a: $(CONFIG)
	@echo "  AR      $@"
	$(Q)[ -d $(@D) ] || mkdir -p $(@D)
	$(Q)$(AR) rcs $@ $(filter %.o,$^)

bin/%.bin: elf/%.elf
	@echo "  OBJCOPY $@"
	$(Q)[ -d $(@D) ] || mkdir -p $(@D)
	$(Q)$(OBJCOPY) -Obinary $< $@

bin/%.hex: elf/%.elf
	@echo "  OBJCOPY $@"
	$(Q)[ -d $(@D) ] || mkdir -p $(@D)
	$(Q)$(OBJCOPY) -Oihex $< $@

obj/%.c.o: %.c $(CONFIG)
	@echo "  CC      $@"
	$(Q)[ -d $(@D) ] || mkdir -p $(@D)
	$(Q)$(CC) -c -o $@ $(CFLAGS) $<

obj/%.c.S: %.c $(CONFIG)
	@echo "  CC      $@"
	$(Q)[ -d $(@D) ] || mkdir -p $(@D)
	$(Q)$(CC) -S -o $@ $(CFLAGS) $<

obj/%.S.o: %.S $(CONFIG)
	@echo "  AS      $@"
	$(Q)[ -d $(@D) ] || mkdir -p $(@D)
	$(Q)$(CC) -c -o $@ $(CFLAGS) $<

obj/%.s.o: %.s $(CONFIG)
	@echo "  AS      $@"
	$(Q)[ -d $(@D) ] || mkdir -p $(@D)
	$(Q)$(CC) -c -o $@ $(CFLAGS) $<

obj/hashprof/%.c.o: %.c $(CONFIG)
	@echo "  CC      $@"
	$(Q)[ -d $(@D) ] || mkdir -p $(@D)
	$(Q)$(CC) -c -o $@ $(CFLAGS) $<

obj/hashprof/%.S.o: %.S $(CONFIG)
	@echo "  AS      $@"
	$(Q)[ -d $(@D) ] || mkdir -p $(@D)
	$(Q)$(CC) -c -o $@ $(CFLAGS) $<

obj/hashprof/%.s.o: %.s $(CONFIG)
	@echo "  AS      $@"
	$(Q)[ -d $(@D) ] || mkdir -p $(@D)
	$(Q)$(CC) -c -o $@ $(CFLAGS) $<

bin-host/%: $(HOST_LIBDEPS) $(CONFIG)
	@echo "  HOST-LD $@"
	$(Q)[ -d $(@D) ] || mkdir -p $(@D)
	$(Q)$(HOST_LD) $(HOST_LDFLAGS) -o $@ $(filter %.o,$^) $(HOST_LDLIBS)

obj-host/%.a: $(CONFIG)
	@echo "  HOST-AR $@"
	$(Q)[ -d $(@D) ] || mkdir -p $(@D)
	$(Q)$(HOST_AR) rcs $@ $(filter %.o,$^)

obj-host/%.c.o: %.c $(CONFIG)
	@echo "  HOST-CC $@"
	$(Q)[ -d $(@D) ] || mkdir -p $(@D)
	$(Q)$(HOST_CC) -c -o $@ $(HOST_CFLAGS) $<

obj-host/%.c.S: %.c $(CONFIG)
	@echo "  HOST-CC $@"
	$(Q)[ -d $(@D) ] || mkdir -p $(@D)
	$(Q)$(HOST_CC) -S -o $@ $(HOST_CFLAGS) $<

obj-host/%.S.o: %.S $(CONFIG)
	@echo "  HOST-AS $@"
	$(Q)[ -d $(@D) ] || mkdir -p $(@D)
	$(Q)$(HOST_CC) -c -o $@ $(HOST_CFLAGS) $<

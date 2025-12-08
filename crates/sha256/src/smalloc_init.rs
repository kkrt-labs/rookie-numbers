use smalloc::Smalloc;
#[global_allocator]
static GLOBAL: Smalloc = Smalloc::new();

#[ctor::ctor]
unsafe fn init_smalloc() {
    GLOBAL.init();
}

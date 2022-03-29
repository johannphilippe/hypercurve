#ifndef RUBY_BRIDGE_H
#define RUBY_BRIDGE_H

#include"core.h"
#include"curve_lib.h"
#include"ruby.h"

VALUE world(VALUE self) {
	printf("HelloWorld\n");
	return Qnil;
}

void Init_ruby_hypercurve() {
	VALUE RHC = rb_define_module("RubyHyperCurve");
	rb_define_singleton_method(RHC, "world", world, 0);
}

#endif // RUBY_BRIDGE_H

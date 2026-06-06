tests = [
	'test_f.f90',
#	'test_p.f90',
#	'test_s.f90',
]

foreach t : tests
  test(
    t,
    executable(
      'test-@0@'.format(t.replace('.f90', '')),
      sources: [t],
      link_with: hompack_lib,
    ),
    suite: 'unit',
    workdir: meson.project_source_root(),
  )
endforeach

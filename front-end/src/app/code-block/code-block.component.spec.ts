import { ComponentFixture, TestBed } from '@angular/core/testing';

import { CodeBlockComponent } from './code-block.component';

describe('CodeBlockComponent', () => {
  let component: CodeBlockComponent;
  let fixture: ComponentFixture<CodeBlockComponent>;

  beforeEach(() => {
    TestBed.configureTestingModule({
      declarations: [CodeBlockComponent]
    });
    fixture = TestBed.createComponent(CodeBlockComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});

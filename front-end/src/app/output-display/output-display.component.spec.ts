import { ComponentFixture, TestBed } from '@angular/core/testing';

import { OutputDisplayComponent } from './output-display.component';
import { MatCardModule } from '@angular/material/card';
import { HttpClientTestingModule } from '@angular/common/http/testing';

describe('OutputDisplayComponent', () => {
  let component: OutputDisplayComponent;
  let fixture: ComponentFixture<OutputDisplayComponent>;

  beforeEach(() => {
    TestBed.configureTestingModule({
      declarations: [OutputDisplayComponent],
      imports: [HttpClientTestingModule,
        MatCardModule],
    });
    fixture = TestBed.createComponent(OutputDisplayComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});

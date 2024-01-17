import { ComponentFixture, TestBed } from '@angular/core/testing';
import { HttpClientTestingModule } from '@angular/common/http/testing';
import { MatCardModule } from '@angular/material/card';
import { MatSnackBarModule } from '@angular/material/snack-bar';

import { OutputDisplayComponent } from './output-display.component';

describe('OutputDisplayComponent', () => {
  let component: OutputDisplayComponent;
  let fixture: ComponentFixture<OutputDisplayComponent>;

  beforeEach(() => {
    TestBed.configureTestingModule({
      declarations: [OutputDisplayComponent],
      imports: [
        HttpClientTestingModule,
        MatCardModule,
        MatSnackBarModule,
      ],
    });
    fixture = TestBed.createComponent(OutputDisplayComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});

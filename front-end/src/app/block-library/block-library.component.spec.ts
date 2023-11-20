import { ComponentFixture, TestBed } from '@angular/core/testing';

import { BlockLibraryComponent } from './block-library.component';
import { MatCardModule } from '@angular/material/card';
import { BlockService } from '../block.service';
import { By } from '@angular/platform-browser';

describe('BlockLibraryComponent', () => {
  let component: BlockLibraryComponent;
  let fixture: ComponentFixture<BlockLibraryComponent>;

  beforeEach(() => {
    TestBed.configureTestingModule({
      declarations: [BlockLibraryComponent],
      imports: [MatCardModule]
    });
    fixture = TestBed.createComponent(BlockLibraryComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });

  describe('AddBlock', () => {

    it ('should add block when button is clicked', () => {
      
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, "addBlock");
      const button = fixture.debugElement.query(By.css("button"));
      button.triggerEventHandler("click", {})
      fixture.detectChanges();
      expect(blockService.addBlock).toHaveBeenCalledOnceWith("LoadData");

    })

  })
});
